
import unittest
import countBarcodeInFastq

class TestCountBarcodeInFastq(unittest.TestCase):

    def setUp(self):
        bs1 = countBarcodeInFastq.BarcodeSet('bs1', 1, 4, 'sample', 56)
        bs1.add('s1', 'NCAT')
        bs1.add('s2', 'TCAT')
        bs1.add('s3','NATT')
        self.barcode_combination = countBarcodeInFastq.BarcodeCombination([bs1])

    def test_count_barcode_in_fastq(self):
        fastq = 'first1M.fastq.gz'
        iden = countBarcodeInFastq.partial(countBarcodeInFastq.identify_barcode_in_read, barcode_combination=self.barcode_combination, direction=0)
        res = countBarcodeInFastq.count_barcode_in_fastq(fastq, 4, iden)
        print(res)


    






