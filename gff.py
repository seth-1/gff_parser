"""
Parser module for gff3 files (dev)
"""

import re
from Bio.SeqFeature import FeatureLocation

STRAND_TRANSLATION = {"+": +1, "-":-1, ".": None}

def strandtranslate(gffStrand):
    """Uses dictionary defined in class to translate gff strand notation to biopython strand notation """
    biopython_strand = STRAND_TRANSLATION[gffStrand]
    return biopython_strand



# TODO CompoundLocation if name exists, match it into compoundlocation from biopython
# TODO generic object for other classes of fields
# TODO compound location for the exons, attributes must equal attributes
# # TODO GENE with compound location,
# TODO throw fasta sequences out....


class GffFeature:
    """main feature object for gff class """
    def __init__(self, location, gff_type, attributes):
        assert isinstance(location, FeatureLocation), type(location)
        self.attributes = attributes
        self.gff_id = self.attributes.pop("ID", None)
        self.location = location
        self.gff_type = str(gff_type)

class Gene(GffFeature):
    """Gene gff feature; Can contain cds and mRNA."""
    def __init__(self, location, gff_type, attributes):
        super().__init__(location, gff_type, attributes)
        self.name = self.attributes.pop("Name", None)
        self.cdss = []
        self.mrnas = []
        self.rrnas = []

    def __repr__(self):
        return "%s, %s" % (self.gff_type, self.location)

    def add_mrna(self, mrna):
        "Adds RNA to parent object"
        assert isinstance(mrna, mRNA)
        self.mrnas.append(mrna)

    def add_rrna(self, rrna):
        "Adds RNA to parent object"
        assert isinstance(rrna, rRNA)
        self.mrnas.append(rrna)

    def add_cds(self, cds):
        "Adds CDS children to the mRNA object"
        assert isinstance(cds, Cds)
        self.cdss.append(cds)


class GffSub(GffFeature):
    """Sub feature class for gff"""
    # NOTE DO I need this class? Wanted to assign parent but parent also needs class check
    def __init__(self, location, gff_type, attributes):
        super().__init__(location, gff_type, attributes)


class mRNA(GffSub):
    """mRNA class, might be split up in later versions"""
    def __init__(self, location, gff_type, parent, attributes):
        super().__init__(location, gff_type, attributes)
        self.name = self.attributes.pop("Name", None)
        self.cdss = []
        self.exons = []

        assert isinstance(parent, Gene), type(parent)
        self.parent = parent

    def __repr__(arg):
        return super().__repr__() + "something something"

    def add_cds(self, cds):
        "Adds CDS children to the mRNA object"
        assert isinstance(cds, Cds)
        self.cdss.append(cds)

    def add_exon(self, exon):
        "Adds CDS children to the mRNA object"
        assert isinstance(exon, Exon), type(exon)
        self.exons.append(exon)


class tRNA(GffSub):
    """tRNA class, might be split up in later versions"""
    def __init__(self, location, gff_type, parent, attributes):
        super().__init__(location, gff_type, attributes)
        self.name = self.attributes.pop("Name", None)
        self.cdss = []
        self.exons = []
        #
        # assert isinstance(parent, Gene), type(parent)
        # self.parent = parent

    def __repr__(arg):
        return super().__repr__() + "something something"

    def add_cds(self, cds):
        "Adds CDS children to the tRNA object"
        assert isinstance(cds, Cds)
        self.cdss.append(cds)

    def add_exon(self, exon):
        "Adds exon to the tRNA object"
        assert isinstance(exon, Exon), type(exon)
        self.exons.append(exon)

class rRNA(GffSub):
    """tRNA class, might be split up in later versions"""
    def __init__(self, location, gff_type, parent, attributes):
        super().__init__(location, gff_type, attributes)
        self.name = self.attributes.pop("Name", None)
        self.cdss = []
        self.exons = []

        assert isinstance(parent, Gene), type(parent)
        self.parent = parent

    def __repr__(arg):
        return super().__repr__() + "something something"

    def add_cds(self, cds):
        "Adds CDS children to the tRNA object"
        assert isinstance(cds, Cds)
        self.cdss.append(cds)

    def add_exon(self, exon):
        "Adds exon to the tRNA object"
        assert isinstance(exon, Exon), type(exon)
        self.exons.append(exon)

class Exon(GffSub):
    """
    Class for exons
    """
    # NOTE: INIT must complete successfully; create

    def __init__(self, location, gff_type, parent, attributes):
        super().__init__(location, gff_type, attributes)
        self.exon_number = self.get_exon_index(self.gff_id)

        assert isinstance(parent, (tRNA, rRNA, mRNA)), "Exon %s must have parent RNA your parent for this exon is %s" % (self.gff_id, type(parent))
        self.parent = parent

    def __repr__(arg):
        return super().__repr__() + "something about exons"

    def get_exon_index(self, arg):
        """extracts exon id """
        exon_reg = re.search("[0-9]+$", arg)
        if not exon_reg:
            raise ValueError("Exons must have integer ID at the end")
            return int(exon_reg.group())


class Cds(GffSub):
    """
    Class for CDS
    """
    def __init__(self, location, gff_type, parent, attributes):
        super().__init__(location, gff_type, attributes)
        assert isinstance(parent,(Gene, mRNA)), "CDS %s must have parent gene or mRNA your parent for this cds with ID %s is %s" % (self.gff_id, parent.gff_id, type(parent))
        self.parent = parent

    def __repr__(arg):
        return super().__repr__() + "something about CDS"


def parse_gff(fpath):
    """
    Parses gff file
    """
    features = {}

    with open(fpath, "r") as target:
        for line in target:
            if "\t" in line:
                delim = "\t"
            else:
                delim = None

            line = line.strip().split(sep=delim, maxsplit=9)

            if line[0].startswith(">"):
                fasta_id = line[0][1:]
                for line in target:
                    if re.search:
                        pass

            if line[0].startswith("#"): # Ignore these lines for now
                continue
            else:
                seqid, source, gff_type, start, end, score, strand, phase, attributes = line

                attributes = {k:v for k, v in[(att.split("=")) for att in attributes.split(";")]}

                location = FeatureLocation(int(start), int(end), strandtranslate(strand))
                parent = attributes.pop("Parent", None)

                if gff_type == "gene":
                    gff_record = Gene(location, gff_type, attributes)

                elif gff_type in ["mRNA", "rRNA", "tRNA", "CDS", "exon"]:
                    parent = features.get(parent, None)

                    if gff_type == "mRNA":
                        gff_record = mRNA(location, gff_type, parent, attributes)
                        features[parent.gff_id].add_mrna(gff_record)

                    if gff_type == "tRNA":
                        gff_record = tRNA(location, gff_type, parent, attributes)
                        # features[parent.gff_id].add_rna(gff_record)

                    if gff_type == "rRNA":
                        gff_record = rRNA(location, gff_type, parent, attributes)
                        features[parent.gff_id].add_rrna(gff_record)

                    if gff_type == "exon":
                        gff_record = Exon(location, gff_type, parent, attributes)
                        features[parent.gff_id].add_exon(gff_record)

                    if gff_type == "CDS":
                        gff_record = Cds(location, gff_type, parent, attributes)
                        features[parent.gff_id].add_cds(gff_record)
                else:
                    continue
                # records.append(gff_record)
                # assert gff_record.gff_id not in features, print(gff_record)
                # TODO DO not understand this assert
                features[gff_record.gff_id] = gff_record

        return features

# gf = parse_gff("GCF_000149205.2_ASM14920v2_genomic.gff")
