from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation
from Bio.Seq import Seq

import copy
import blastn
import sequences


class AlignedFeature(SeqFeature):
    def __init__(self, hsp, host, donor):
        self.hsp = hsp
        self.donor = donor
        self.donor_location = blastn.hsp_location(hsp, target=True)
        location = blastn.hsp_location(hsp, target=False)
        qualifiers = {
            "label": self.donor.name,
            "donor_start": self.donor_location.start,
            "donor_end": self.donor_location.end,
            "donor_strand": self.donor_location.strand,
        }
        super().__init__(location, "misc_feature", self.donor.id, qualifiers)
        

class InsertedFeature(AlignedFeature):
    def __init__(self, *args):
        super().__init__(*args)

def insertion_search(reads, donors):
    for alignment in blastn.blastn(reads, donors, word_size=12):
        for hit in alignment:
            for hsp in hit:
                read = next(x for x in reads if x.id == hsp.query.id)
                donor = next(x for x in donors if x.id == hsp.target.id)
                read.features.append(InsertedFeature(hsp, read, donor))

class GenomeFeature(AlignedFeature):
    def __init__(self, *args):
        super().__init__(*args)

def genome_search(reads, donors):
    masked_reads = [read.mask_insertion() for read in reads]
    for alignment in blastn.blastn(masked_reads, donors):
        for hit in alignment:
            for hsp in hit:
                read = next(x for x in reads if x.id == hsp.query.id)
                donor = next(x for x in donors if x.name == hsp.target.name)
                print(f"Genome HSP in {read.name}")
                read.features.append(GenomeFeature(hsp, read, donor))

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
    def aligned_features(self):
        "Features of this read that come from alignments."
        return [f for f in self.features if isinstance(f, AlignedFeature)]

    @property
    def inserted_features(self):
        "Features of the read that come from transposon insertions."
        return [f for f in self.features if isinstance(f, InsertedFeature)]

    @property
    def has_insertion(self):
        "Returns True if this read has alignments from a transposon insertion."
        return bool(self.inserted_features)
    
    @property
    def hsps(self):
        "All the HSPs associated with this read."
        return [f.hsp for f in self.aligned_features]

    def aligned_features_from(self, donor):
        "Features of this read that come from alignments with donor."
        return [f for f in self.aligned_features if f.donor.id == donor.id]

    def hsps_from(self, donor):
        "The HSPs associated with this read that come from donor."
        return [f.hsp for f in self.aligned_features_from(donor)]

    def has_hsp_from(self, donor):
        "Returns True if this read has HSPs from donor"
        return bool(self.hsps_from(donor))

    def transposon_search(self, transposon_sequences):
        return transposon_search([self], transposon_sequences)

    def mask_insertion(self):
        read = copy.deepcopy(self)
        if self.has_insertion:
            mask = sum(i.location for i in self.inserted_features)
            seq = ("N" if i in mask else a for (i, a) in enumerate(self))
            seq = Seq("".join(seq))
            read = copy.deepcopy(self)
            read.seq = seq
        return read

    def merge_features_with(self, other):
        for feature in other.features:
            if feature not in self.features:
                self.features.append(feature)
            
        
