from Bio.SeqRecord import SeqRecord

import blastn

def transposon_search(insertions, transposons):
    blastn.align(insertions, transposons)
    for insertion in insertions:
        transposon_features = []
        for transposon in transposons:
            transposon_features += insertion.hsp_features_from(transposon)

        if transposon_features:
            scorefun = lambda f: blastn.hsp_score(f.hsp)
            best = sorted(transposon_features, key=scorefun)[-1]
            new_features = [f for f in insertion.features if f not in transposon_features]
            insertion.features = new_features
            insertion.features.append(best)
            insertion.transposon = best.hsp.target

    return None
        

class Insertion(SeqRecord):
    def __init__(self, seqrecord, **kwargs):
        kwargs = kwargs | {
            "id": seqrecord.id,
            "name": seqrecord.name,
            "description": seqrecord.description,
            "dbxrefs": seqrecord.dbxrefs,
            "features": seqrecord.features,
            "annotations": seqrecord.annotations,
            "letter_annotations": seqrecord.letter_annotations,
        }
        super().__init__(seqrecord.seq, **kwargs)
        self.transposon = None
        self.position = None

    @property
    def hsp_features(self):
        "Features of this insertion which are associated with HSPs"
        return [f for f in self.features if hasattr(f, "hsp")]

    def hsp_features_from(self, target):
        "Features of this insertion which are associated with HSPs paired with target"
        return [f for f in self.hsp_features if f.hsp.target.name == target.name]

    @property
    def hsps(self):
        "The HSPs associated with this insertion"
        return [f.hsp for f in self.hsp_features]

    def hsps_from(self, target):
        "The HSPs associated with this insertion paired with target"
        return [h for h in self.hsps if h.target.name == target.name]

    def has_hsp_from(self, target):
        "Returns True if this insertion has HSPs paired with target"
        return bool(self.hsps_from(target))

    @property
    def has_transposon(self):
        "Returns True if this insertion is known to include transposon sequence"
        return self.transposon is not None

    def transposon_search(self, transposon_sequences):
        return transposon_search([self], transposon_sequences)
        


            
        
