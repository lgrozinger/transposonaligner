from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation
from Bio.Seq import Seq

import blastn
import sequences

class AlignedFeature(SeqFeature):
    pass

class TransposonFeature(AlignedFeature):
    def __init__(self, hsp, parent):
        query_location = blastn.hsp_to_location(hsp, target=False)
        target_location = blastn.hsp_to_location(hsp, target=True)
        self.hsp = hsp

class InsertedFeature(AlignedFeature):
    def __init__(self, hsp, host, donor):
        self.hsp = hsp
        self.donor = donor
        self.donor_location = blastn.hsp_location(hsp, target=True)
        location = blastn.hsp_location(hsp, target=False)
        qualifiers = {"label": self.donor.name}
        super().__init__(location, "misc_feature", self.donor.id, qualifiers)

def transposon_search(reads, donors):
    for alignment in blastn.blastn(reads, donors, word_size=12):
        for hit in alignment:
            for hsp in hit:
                read = next(x for x in reads if x.id == hsp.query.id)
                donor = next(x for x in donors if x.id == hsp.target.id)
                read.features.append(InsertedFeature(hsp, read, donor))

class Read(SeqRecord):
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
    def insertion_features(self):
        return [f for f in self.features if isinstance(f, InsertionFeature)]
    
    @property
    def has_transposon(self):
        "Returns True if this insertion is known to include transposon sequence"
        return self.transposon is not None

    def transposon_search(self, transposon_sequences):
        return transposon_search([self], transposon_sequences)

    def mask_insertion(self):
        mask = sum(i.location for i in self.insertion_features)
        seq = ("N" if i in mask else a for (i, a) in enumerate(self))
        seq = Seq("".join(seq))
        read = copy.deepcopy(self)
        read.seq = seq
        return read

    def between(self, x, y):
        record = sequences.region_between(self, x, y)
        return Read(record)

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
            
        
