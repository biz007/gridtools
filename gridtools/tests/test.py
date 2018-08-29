from unittest import TestCase

from gridtools import grid, main

class Test(TestCase):
    def test_model(self):
        s = main.evalModel
        print(s)