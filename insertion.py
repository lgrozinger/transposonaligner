from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation
from Bio.Seq import Seq

import blastn
import sequences

def query(hit, insertions):
    return next(filter(lambda x: x.id == hit[0].query.id, insertions))

def target(hit, insertions):
    return next(filter(lambda x: x.name == hit.target.name, insertions))

def transposon_search(insertions, transposons):
    for alignment in blastn.blastn(insertions, transposons, word_size=8):
        for hit in alignment:
            insertion = query(hit, insertions)
            transposon = target(hit, transposons)
            insertion.features.append(TransposonFeature(hit, insertion))

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

class AlignedFeature(SeqFeature):
    pass

class TransposonFeature(AlignedFeature):
    def __init__(self, hsp, parent):
        query_location = blastn.hsp_to_location(hsp, target=False)
        target_location = blastn.hsp_to_location(hsp, target=True)
        self.hsp = hsp

        import pdb; pdb.set_trace()
                
class Insertion(SeqRecord):
    def __init__(
            self,
            seqrecord_or_seq,
            id = "<unknown id>",
            name = "<unknown name>",
            description = "<unknown description>",
            dbxrefs = None,
            features = None,
            annotations = None,
            letter_annotations = None,
            ):

        if isinstance(seqrecord_or_seq, SeqRecord):
            seq = seqrecord_or_seq.seq
            id = seqrecord_or_seq.id
            name = seqrecord_or_seq.name
            description = seqrecord_or_seq.description
            dbxrefs = seqrecord_or_seq.dbxrefs
            features = seqrecord_or_seq.features
            annotations = seqrecord_or_seq.annotations
            letter_annotations = seqrecord_or_seq.letter_annotations
        elif isinstance(seqrecord_or_seq, Seq):
            seq = seqrecord_or_seq
        else:
            raise TypeError(f"Cannot construct Insertion from {type(seqrecord_or_seq)}")
            
        super().__init__(
            seq,
            id,
            name,
            description,
            dbxrefs,
            features,
            annotations,
            letter_annotations
        )
        
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

    def between(self, x, y):
        record = sequences.region_between(self, x, y)
        return Insertion(record)

    @property
    def last_aligned_base(self):
        last = 0
        for feature in self.features:
            last = max(last, feature.location.start, feature.location.end)
        return last

    @property
    def suffix(self):
        return self.between(self.last_aligned_base, len(self) + 1)

    def merge_features_with(self, other):
        for feature in other.features:
            if feature not in self.features:
                self.features.append(feature)
            
        
