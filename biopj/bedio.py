from biopj import genomicelements


class BedFile:
    """
    Class to read and write BED files
    """

    def __init__(self, file_something):
        """
        Either a string (filename) or file object has to be provided as file_something.
        If filename is an existing file, it will be opened for reading otherwise the file will be created and opened for
        writing.

        :type file_something: str/File
        """
        if type(file_something) is str:
            self.filename = file_something
            try:
                self.file = open(file_something, 'r+')
                self.mode_reading = True
            except IOError as e:
                if e.errno == 2:
                    self.file = open(file_something, 'w+')
                    self.mode_reading = False
                else:
                    raise e

        elif hasattr(file_something, 'read'):
            self.file = file_something
            self.mode_reading = True
        else:
            raise Exception("No filename or file object provided")

    def __iter__(self):
        if self.mode_reading:
            for line in self.file:
                yield BedLine(line)
        else:
            raise Exception("File not opened for reading")


    def write(self, something):
        """
        Write a line to the BED file. Parameter something will be converted to string.
        If something is a genomicelements.Transcript object, it will be converted to
        a BedLine object

        :param something: object
        """
        if type(something) is genomicelements.Transcript:
            bl = BedLine(something)
            self.file.write(str(bl) + '\n')
        else:
            self.file.write(str(something) + '\n')

    def close(self):
        self.file.close()


class BedLine:
    """
    Line in a BED file. Can be converted to and from genomicelements.Transcript objects
    :type blocks: list[BedBlock]
    """
    
    blocks = []

    def __init__(self, info: object):
        """
        :param info: str or Transcript
        """
        if type(info) is genomicelements.Transcript:
            self.fromTranscript(info)
        elif type(info) is str:
            self.fromString(info)
        else:
            raise Exception("Cannot initiate BedFile object with type " + type(info))

    def fromTranscript(self, transcript):
        """
        :type transcript: genomicelements.Transcript
        :param transcript: Transcript
        """
        self.chrom = transcript.chrom
        self.chrom_start = transcript.start - 1
        self.chrom_end = transcript.end
        self.strand = transcript.strand
        self.name = transcript.name

        # if these fields are absent, BedLine will not be complete
        if len(transcript.exons) > 0:
            self.score       = 0
            self.thick_start = self.chrom_start
            self.thick_end   = self.chrom_end
            self.item_rgb    = '0,0,0'

        for exon in transcript.exons:
            block = BedBlock(transcript.chrom, exon.start + 1, exon.end, transcript.strand)
            self.blocks.append(block)

    def fromString(self, line):
        """
        :param line: str
        """
        l = line.rstrip().split('\t')

        if len(l) < 3:
            raise Exception("Improperly formatted BED line '%s'"%line)

        self.blocks = []

        self.chrom = l[0]
        self.chrom_start = int(l[1])  # 0-based start
        self.chrom_end = int(l[2])  # 1-based end

        if len(l) > 3:
            self.name = l[3]

            if len(l) > 4:
                self.score = int(l[4])

                if len(l) > 5:
                    self.strand = l[5]

                    if len(l) > 6:
                        self.thick_start = int(l[6])
                        self.thick_end = int(l[7])

                        if len(l) > 8:
                            self.item_rgb = l[8]

                            if len(l) > 9:
                                block_count = int(l[9])
                                block_sizes = l[10].split(',')
                                block_starts = l[11].split(',')

                                for i in range(0, block_count):
                                    block_start = int(block_starts[i]) + self.chrom_start
                                    block_end = block_start + int(block_sizes[i])
                                    block = BedBlock(self.chrom, block_start, block_end, self.strand)

                                    self.blocks.append(block)

    def __str__(self):
        r_val = '%s\t%s\t%s' % (self.chrom, self.chrom_start, self.chrom_end)

        if hasattr(self, 'name'):
            r_val += '\t%s' % (self.name)

            if hasattr(self, 'score'):
                r_val += '\t%s' % (self.score)

                if hasattr(self, 'strand'):
                    r_val += '\t%s' % (self.strand)

                    if hasattr(self, 'thick_start'):
                        r_val += '\t%s\t%s' % (self.thick_start, self.thick_end)

                        if hasattr(self, 'item_rgb'):
                            r_val += '\t%s' % (self.item_rgb)

                            if len(self.blocks) > 0:
                                block_count = len(self.blocks)
                                block_sizes = ','.join(map((lambda x: str(x.end - x.start)), self.blocks)) + ','
                                block_starts = ','.join(
                                    map((lambda x: str(x.start - self.chrom_start)), self.blocks)) + ','
                                r_val += '\t%s\t%s\t%s' % (block_count, block_sizes, block_starts)

        return r_val

    def toGenomicElement(self):
        return genomicelements.GenomicElement(self.chrom, self.chrom_start - 1, self.chrom_end, self.strand)

    def toTranscript(self):
        t = genomicelements.Transcript(self.chrom, self.chrom_start - 1, self.chrom_end, self.strand, self.name)
        
        for b in self.blocks:
            e = genomicelements.Exon(b.chrom, b.start - 1, b.end, b.strand, self.name)
            t.exons.append(e)

        return t


class BedBlock:
    def __init__(self, chrom, start, end, strand):
        """

        """
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand

    def toGenomicElement(self):
        return genomicelements.GenomicElement(self.chrom, self.start - 1, self.end, self.strand)
