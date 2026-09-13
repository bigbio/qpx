"""Order-independent lookup of protein groups already inferred by OpenMS."""

from __future__ import annotations

from dataclasses import dataclass, field


def identification_identifier(identification) -> str:
    """Read the source key from a streamed record or a pyOpenMS identification."""
    identifier = getattr(identification, "identifier", None)
    return identification.getIdentifier() if identifier is None else identifier


@dataclass
class ProteinGroupIndex:
    """Keep full memberships and every accession's possible groups."""

    by_membership: dict[frozenset[str], tuple[str, ...]] = field(default_factory=dict)
    by_accession: dict[str, set[frozenset[str]]] = field(default_factory=dict)
    by_identification: dict[str, ProteinGroupIndex] = field(default_factory=dict)

    @classmethod
    def from_groups(cls, groups) -> ProteinGroupIndex:
        """Index memberships while preserving each producer-supplied leader."""
        index = cls()
        for group in groups:
            key = frozenset(group)
            if not key:
                continue
            index.by_membership.setdefault(key, tuple(group))
            for accession in key:
                index.by_accession.setdefault(accession, set()).add(key)
        return index

    def resolve(self, accessions, identifier=None) -> tuple[str, ...] | None:
        """Match the complete group first, else require one linked group.

        Evidence may include proteins excluded by upstream inference, so do not
        require the inferred group to contain every possible sequence match.
        A known identification's groups take precedence over the merged index.
        """
        index = self.by_identification.get(identifier) if identifier and self.by_identification else self
        if index is None or not accessions:
            return None
        key = frozenset(accessions)
        if key in index.by_membership:
            return self.by_membership[key]
        candidates = set().union(*(index.by_accession.get(acc, ()) for acc in key))
        if len(candidates) == 1:
            return self.by_membership[next(iter(candidates))]
        return None

    def unambiguous_accessions(self) -> dict[str, list[str]]:
        """Legacy accession lookup, excluding proteins shared across groups."""
        return {
            acc: list(self.by_membership[next(iter(groups))]) for acc, groups in self.by_accession.items() if len(groups) == 1
        }
