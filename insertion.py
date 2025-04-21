from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation
from Bio.Seq import Seq

import copy
import pandas

from blastn import blastn
from alignedfeature import InsertedFeature
from alignedfeature import GenomeFeature


def insertion_search(reads, donors):
    options = {"evalue": 0.01}
    for alignment in blastn(reads, donors, **options):
        for hit in alignment:
            for hsp in hit:
                read = next(x for x in reads if x.id == hsp.query.id)
                donor = next(x for x in donors if x.id == hsp.target.id)
                insert = InsertedFeature(hsp, read, donor)
                read.features.append(insert)
                read.features += insert.donor_features

def genome_search(reads, donors):
    masked_reads = [read.mask_insertion() for read in reads]
    options = {"evalue": 0.01, "word_size": 10}
    for alignment in blastn(masked_reads, donors, **options):
        for hit in alignment:
            for hsp in hit:
                read = next(x for x in reads if x.id == hsp.query.id)
                donor = next(x for x in donors if x.name == hsp.target.name)
                genome = GenomeFeature(hsp, read, donor)
                read.features.append(genome)
                read.features += genome.donor_features

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
    def transposon(self):
        "Returns the supposed transposon system inserted into this read."
        if len(self.inserted_features) > 1:
            raise TypeError("Must run choose_insertion first")
        return self.inserted_features[0] if self.has_insertion else None

    @property
    def genome_features(self):
        "Features of the read that come from a genome."
        return [f for f in self.features if isinstance(f, GenomeFeature)]

    @property
    def has_insertion(self):
        "Returns True if this read has alignments from a transposon insertion."
        return bool(self.inserted_features)

    @property
    def has_genome(self):
        "Returns True if this read has alignments from a genome."
        return bool(self.genome_features)

    @property
    def hsps(self):
        "All the HSPs associated with this read."
        return [f.hsp for f in self.aligned_features]

    @property
    def default_datarow(self):
        return {
            "name": self.id,
            "read length": len(self),
            
            "transposon": None,
            "transposon evalue": None,
            "transposon bit score": None,
            "transposon identity": None,
            
            "genome": None,
            "genome evalue": None,
            "genome bit score": None,
            "genome identity": None,
            "genome length": None,
            "genome mismatch": None,
            "genome gaps": None,
            "genome start": None,
            "genome end": None,
            "genome strand": None,
            
            "insertion": None,
            "insertion name": None,
            "insertion type": None,
            "insertion strand": None,
            "insertion product": None,
    }
    
    @property
    def dataframe(self):
        rows = []
        common = self.default_datarow
        if self.has_insertion:
            common = common | {
                "transposon": self.transposon.donor.name,
                "transposon evalue": self.transposon.evalue,
                "transposon bit score": self.transposon.bitscore,
                "transposon identity": self.transposon.identity,
                "transposon length": int(self.transposon.length),
                "transposon start": int(self.transposon.location.start),
                "transposon end": int(self.transposon.location.end),
            }
            
        if self.has_genome:
            for candidate in self.genome_features:
                rows.append(common | candidate.dataframe)
        else:
            rows.append(common)

        return pandas.DataFrame(rows)

    def aligned_features_from(self, donor):
        "Features of this read that come from alignments with donor."
        return [f for f in self.aligned_features if f.donor.id == donor.id]

    def hsps_from(self, donor):
        "The HSPs associated with this read that come from donor."
        return [f.hsp for f in self.aligned_features_from(donor)]

    def has_hsp_from(self, donor):
        "Returns True if this read has HSPs from donor"
        return bool(self.hsps_from(donor))

    def insertion_search(self, inserted_sequences):
        return insertion_search([self], inserted_sequences)

    def mask_insertion(self):
        read = copy.deepcopy(self)
        if self.has_insertion:
            mask = sum(i.location for i in self.inserted_features)
            seq = ("N" if i in mask else a for (i, a) in enumerate(self))
            seq = Seq("".join(seq))
            read = copy.deepcopy(self)
            read.seq = seq
        return read

    def choose_insertion(self):
        if self.has_insertion:
            best = sorted(self.inserted_features, key=lambda x: len(x))[-1]
            for insertion in self.inserted_features:
                self.features.remove(insertion)
                for feature in insertion.donor_features:
                    self.features.remove(feature)
            self.features.append(best)
            self.features += best.donor_features

    def choose_genome(self):
        if self.has_genome:
            scorer = lambda x: (min(x.location.start, x.location.end), -len(x))
            best_score = min(scorer(x) for x in self.genome_features)
            best = [x for x in self.genome_features if scorer(x) == best_score]
            for genome in self.genome_features:
                self.features.remove(genome)
                for feature in genome.donor_features:
                    self.features.remove(feature)
            for genome in best:
                self.features.append(genome)
                self.features += genome.donor_features

