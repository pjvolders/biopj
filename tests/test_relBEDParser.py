from unittest import TestCase

from biopj import relbed
from biopj.bedio import BedFile

class TestRelBEDParser(TestCase):

    def test___init__(self):
        rb = relbed.RelBEDParser('assets/test_file.bed')

        self.assertEqual(rb.index['NNT-AS1:9'].chrom, 'chr5')
        self.assertEqual(rb.index['PITPNM2-AS1:2'].strand, '+')

    # def test_parse_line(self):
    #     self.fail()

    def test_parse_line(self):
        rb = relbed.RelBEDParser('assets/test_file.bed')
        rbf = BedFile('assets/test_relfile.bed')
        rel_items = []
        for i in rbf:
            rel_items.append(i)
        rbf.close()

        # single exon +
        parsed_line = rb.parse_line(rel_items[0])
        self.assertEqual('chr12', parsed_line.chrom)
        self.assertEqual(123081385, parsed_line.chrom_start)
        self.assertEqual(123081483, parsed_line.chrom_end)

        # single exon -
        #PIK3CD-AS2:1	5	634	feature_2
        #chr1	9679799	9687574	PIK3CD-AS2:1	0	-	9679799	9679799	0,0,0	1	7775,	0,
        parsed_line_2 = rb.parse_line(rel_items[1])
        self.assertEqual('chr1', parsed_line_2.chrom)
        self.assertEqual(9686940, parsed_line_2.chrom_start)
        self.assertEqual(9687569, parsed_line_2.chrom_end)

        # multi exon +
        #DUBR:2	10	1000	feature_3
        #chr3	107240691	107248348	DUBR:2	0	+	107240691	107240691	0,0,0	2	396,2341,	0,5316,
        parsed_line_3 = rb.parse_line(rel_items[2])
        self.assertEqual('chr3', parsed_line_3.chrom)
        self.assertEqual(107240701, parsed_line_3.chrom_start)
        self.assertEqual(107246611, parsed_line_3.chrom_end)

        self.assertEqual(107241087, parsed_line_3.blocks[0].end)
        self.assertEqual(107246007, parsed_line_3.blocks[1].start)
        self.assertEqual(107246611, parsed_line_3.blocks[1].end)

        # muti exon -
        # NNT-AS1:9	10	300	feature_4
        # chr5	43578336	43602937	NNT-AS1:9	0	-	43578336	43578336	0,0,0	3	243,100,160,	0,5130,24441,
        parsed_line_4 = rb.parse_line(rel_items[3])
        self.assertEqual('chr5', parsed_line_4.chrom)
        self.assertEqual(43578539, parsed_line_4.chrom_start)
        self.assertEqual(43602927, parsed_line_4.chrom_end)

        self.assertEqual(43578579, parsed_line_4.blocks[0].end)
        self.assertEqual(43602777, parsed_line_4.blocks[2].start)
        self.assertEqual(43602927, parsed_line_4.blocks[2].end)
