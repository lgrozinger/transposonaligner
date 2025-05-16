from Bio.SeqFeature import SeqFeature

from tnatlas.blastn import hsp_location
from tnatlas.annotate import translated_features_in


class AlignedFeature(SeqFeature):
    def __init__(self, hsp, host, donor, **kwargs):
        location = hsp_location(hsp, target=False)
        super().__init__(location, "misc_feature", donor.id)
        
        self.hsp = hsp
        self.donor = donor
        self.donor_location = hsp_location(hsp, target=True)
        self.qualifiers = {
            "label": self.donor.name,
            "donor_start": self.donor_location.start,
            "donor_end": self.donor_location.end,
            "donor_strand": "plus" if self.donor_location.strand > 0 else "minus",
            "evalue": self.hsp.annotations["evalue"],
        }
        self.donor_features = translated_features_in(
            self.donor_location,
            self.location,
            self.donor,
        )

        self.bitscore = self.hsp.annotations["bit score"]
        self.evalue = self.hsp.annotations["evalue"]
        self.positive = self.hsp.annotations["positive"]
        self.length = self.hsp.length
        self.identity = self.hsp.annotations["identity"] / self.length
        self.mismatch = self.length - self.hsp.annotations["identity"]
        self.gaps = self.hsp.annotations["gaps"]
        self.prefix_length = kwargs.get("prefix_length", None)

    def sequence_prefix(self, n=None):
        if n is None:
            n = self.prefix_length            
        start = min(self.hsp.coordinates[1])
        finish = min(max(self.hsp.coordinates[1]), start + n)
        return self.hsp.query.seq[start:finish]

    def alignment_prefix(self, n=None):
        if n is None:
            n = self.prefix_length
        midlines = self.hsp.annotations["midline"]
        finish = min(len(midlines), n)
        return midlines[:finish].replace(" ", ".")
        
class InsertedFeature(AlignedFeature):
    def __init__(self, *args, **kwargs):
        super().__init__(*args)

class GenomeFeature(AlignedFeature):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def loci(self):
        loci = [f.qualifiers.get("locus_tag", []) for f in self.donor_features]
        loci = [l for ls in loci for l in ls]
        return list(set(loci))

    @property
    def interesting_features(self):
        interesting = lambda f: f.type not in ["source", "gene"]
        return [f for f in self.donor_features if interesting(f)]
    
    @property
    def disrupted(self):
        features = self.interesting_features
        startof = lambda x: min(x.location.start, x.location.end)
        return [f for f in features if startof(f) == self.location.start]

    @property
    def disrupted_dataframe(self):
        features = self.disrupted
        loci = [f.qualifiers.get("locus_tag", [None]) for f in features]
        loci = [l for ls in loci for l in ls]
        genes = [f.qualifiers.get("gene", [None]) for f in features]
        genes = [g for gs in genes for g in gs]
        types = [f.type for f in features]
        strands = [f.location.strand for f in features]
        strands = ["plus" if s > 0 else "minus" for s in strands]
        products = [f.qualifiers.get("product", [None]) for f in features]
        products = [p for ps in products for p in ps]
        return {
            "insertion locus": loci if loci else None,
            "insertion gene": genes if genes else None,
            "sequence type": types if features else None,
            "insertion strand": strands if strands else None,
            "insertion product": products if products else None,
        }
    
    @property
    def dataframe(self):
        return self.disrupted_dataframe | {
            "genome": self.donor.id,
            "genome alignment evalue": self.evalue,
            "genome alignmenr bit score": self.bitscore,
            "genome alignment identity": self.identity,
            "genome alignment length": int(self.length),
            "genome alignment mismatch": self.mismatch,
            "genome alignment gaps": self.gaps,
            "genome start": int(self.donor_location.start),
            "genome end": int(self.donor_location.end),
            "genome strand": "plus" if self.donor_location.strand > 0 else "minus",
            "genome loci": self.loci,
            "genome prefix": self.sequence_prefix(),
            "genome alignment": self.alignment_prefix(),
        }
