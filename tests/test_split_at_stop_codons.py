from unittest import TestCase
from biopj.proteomics.split_at_stop_codons import find_splits, split_sequence


class Test(TestCase):
    def test_find_splits(self):
        self.assertListEqual(find_splits("HELLO*WORLD"), [0, 6, 11])
        self.assertListEqual(find_splits("*HELLO*WORLD*"), [0, 1, 7, 13])

    def test_split_sequence(self):
        splits = list(split_sequence("HELLO*WORLD"))
        self.assertListEqual(splits, [("HELLO", 0, 6), ("WORLD", 6, 11)])

        splits = list(split_sequence("*HELLO*WORLD*"))
        self.assertListEqual(splits, [("HELLO", 1, 7), ("WORLD", 7, 13)])
        #self.fail()
