from unittest import TestCase

from biopj.bedio import BedFile

class TestBedFile(TestCase):
    def test___iter__(self):
        bf = BedFile('assets/test_file.bed')
        n = 0
        for bl in bf:
            n += 1
        bf.close()
        self.assertEqual(n, 9, 'BED file has 9 lines')
        self.assertEqual(bl.name, 'DUBR:2', 'Last line in BED file is DUBR:2')
