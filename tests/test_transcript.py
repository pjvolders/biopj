from unittest import TestCase

from biopj.genomicelements import Transcript, Exon


class TestTranscript(TestCase):
    def test_size_none(self):
        a = Transcript("chr1", 83801516, 83860546, "-", "LINC01725")
        self.assertEqual(a.size, None)

    def test_size(self):
        a = Transcript("chr1", 83801516, 83860546, "-", "LINC01725")
        a1 = Exon("chr1", 83801516, 83803251, "+", a)
        a.exons.append(a1)
        a2 = Exon("chr1", 83849907, 83850022, "+", a)
        a.exons.append(a2)
        a3 = Exon("chr1", 83860408, 83860546, "+", a)
        a.exons.append(a3)

        self.assertEqual(a.size, 1991)
