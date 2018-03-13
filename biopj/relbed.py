#!/usr/bin/env python

# import modules used here
from biopj import bedio
from biopj.genomicelements import Transcript


class RelBEDParser(object):
    """docstring for RelBEDParser
    :type index: dict[str,bedio.BedLine]
    """
    def __init__(self, bed_file):
        super(RelBEDParser, self).__init__()
        self.bed_file = bedio.BedFile(bed_file)
        self.index = {}

        for bl in self.bed_file:
            if bl.name in self.index:
                raise Exception("Name is not unique in BED file")
            else:
                self.index[bl.name] = bl
        self.bed_file.close()

    def parse_line(self, relbed_line):
        """
        Convert a BedLine object to the coordinate system in the bed_file
        :type relbed_line: bedio.BedLine
        :param relbed_line:
        :return: bedio.BedBLock
        """
        if relbed_line.chrom not in self.index:
            raise Exception("Parent transcript does not exist in index")
        else:
            bed_line = self.index[relbed_line.chrom]

            if bed_line.strand == '+':
                if len(bed_line.blocks) <= 1:
                    abs_start = bed_line.chrom_start + relbed_line.chrom_start
                    abs_end = bed_line.chrom_start + relbed_line.chrom_end

                else:
                    for block in bed_line.blocks:
                        exon_size = block.end-block.start
                        if relbed_line.chrom_start >= exon_size:
                            relbed_line.chrom_start -= exon_size
                        else:
                            abs_start = block.start + relbed_line.chrom_start
                            break

                    for block in bed_line.blocks:
                        exon_size = block.end-block.start
                        if relbed_line.chrom_end > exon_size:
                            relbed_line.chrom_end -= exon_size
                        else:
                            abs_end = block.start + relbed_line.chrom_end
                            break
            else:
                if len(bed_line.blocks) <= 1:
                    abs_start = bed_line.chrom_end - relbed_line.chrom_end
                    abs_end = bed_line.chrom_end - relbed_line.chrom_start

                else:
                    for block in reversed(bed_line.blocks):
                        exon_size = block.end-block.start
                        if relbed_line.chrom_start >= exon_size:
                            relbed_line.chrom_start -= exon_size
                        else:
                            abs_end = block.end - relbed_line.chrom_start
                            break

                    for block in reversed(bed_line.blocks):
                        exon_size = block.end-block.start
                        if relbed_line.chrom_end > exon_size:
                            relbed_line.chrom_end -= exon_size
                        else:
                            abs_start = block.end - relbed_line.chrom_end
                            break

            # get the blocks
            abs_blocks = list(filter(lambda b: b.end >= abs_start and b.start <= abs_end, bed_line.blocks))
            if len(abs_blocks) > 0:
                abs_blocks[0].start = abs_start
                abs_blocks[-1].end = abs_end

            tr = Transcript(bed_line.chrom, abs_start+1, abs_end, bed_line.strand, relbed_line.name)
            bl_abs = bedio.BedLine(tr)
            bl_abs.blocks = abs_blocks
            return bl_abs



