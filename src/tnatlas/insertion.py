from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation
from Bio.Seq import Seq

import copy
import pandas

from tnatlas.blastn import blastn
from tnatlas.alignedfeature import InsertedFeature
from tnatlas.alignedfeature import GenomeFeature


def insertion_search(reads, donors, **kwargs):
    for alignment in blastn(reads, donors, **kwargs):
        for hit in alignment:
            for hsp in hit:
                read = next(x for x in reads if x.id == hsp.query.id)
                donor = next(x for x in donors if x.id == hsp.target.id)
                insert = InsertedFeature(hsp, read, donor)
                read.features.append(insert)
                read.features += insert.donor_features

def genome_search(reads, donors, prefix_length, **kwargs):
    masked_reads = [read.mask_insertion() for read in reads]
    for alignment in blastn(masked_reads, donors, **kwargs):
        for hit in alignment:
            for hsp in hit:
                read = next(x for x in reads if x.id == hsp.query.id)
                donor = next(x for x in donors if x.name == hsp.target.name)
                genome = GenomeFeature(hsp, read, donor, prefix_length=prefix_length)
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
        return {"name": self.id, "read length": len(self)}

    @property
    def genome_offset(self):
        offset = None
        if self.has_genome and self.has_insertion:
            starts = [min(x.location.start, x.location.end) for x in self.genome_features]
            if len(set(starts)) > 1:
                raise TypeError("Must run choose_genome before asking for offset.")

            location = self.genome_features[0].location
            offset = min(location.start, location.end) - self.transposon.location.end
        return offset

    @property
    def gap(self):
        if self.has_genome and self.has_insertion:
            x = self.transposon.location.end
            y = min(genome.location.start for genome in self.genome_features)
            return self[x:y]
        else:
            return None
    
    @property
    def dataframe(self):
        rows = []
        common = {
            "name": self.id,
            "read length": len(self),
            "multiple candidates": len(self.genome_features) > 1,
        }
        if self.has_insertion:
            common["transposon"] = self.transposon.donor.name
            common["transposon evalue"] = self.transposon.evalue
            common["transposon bit score"] = self.transposon.bitscore
            common["transposon identity"] = self.transposon.identity
            common["transposon length"] = int(self.transposon.length)
            common["transposon start"] = int(self.transposon.location.start)
            common["transposon end"] = int(self.transposon.location.end)
            common["genome offset"] = self.genome_offset
            common["gap sequence"] = self.gap.seq if self.genome_offset is not None else None
            
        if self.has_genome:
            for candidate in self.genome_features:
                data = candidate.dataframe
                rows.append(common | data)
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
        "Keep only the best insertion vector alignment, remove the rest."
        if self.has_insertion:
            def sortkey(x):
                return (len(x), x.location.start, -x.location.end)
            best = sorted(self.inserted_features, key=sortkey)[-1]
            for insertion in self.inserted_features:
                self.features.remove(insertion)
                for feature in insertion.donor_features:
                    self.features.remove(feature)
            self.features.append(best)
            self.features += best.donor_features

    def choose_genome(self):
        "Remove all but the first (leftmost) genome alignments."
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

