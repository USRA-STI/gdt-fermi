# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Contract No.: CA 80MSFC17M0022
# Contractor Name: Universities Space Research Association
# Contractor Address: 7178 Columbia Gateway Drive, Columbia, MD 21046
#
# Copyright 2017-2022 by Universities Space Research Association (USRA). All rights reserved.
#
# Developed by: William Cleveland and Adam Goldstein
#               Universities Space Research Association
#               Science and Technology Institute
#               https://sti.usra.edu
#
# Developed by: Daniel Kocevski
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing permissions and limitations under the
# License.
#

import os
import numpy as np
import unittest

from gdt.missions.fermi.gbm.collection import *

class MyDataObject():
    def __init__(self, detector, value):
        self._detector = detector
        self._value = value
    
    @property
    def detector(self):
        return self._detector
    
    @property
    def value(self):
        return self._value
    
    def set_value(self, new_val):
        self._value = new_val

    def value_as_str(self):
        return str(self._value)   
    
    def add_int(self, val=1):
        return self._value + val
    

class TestGbmDetectorCollection(unittest.TestCase):
    
    def setUp(self):
        self.b0 = MyDataObject('b0', 0) 
        self.b1 = MyDataObject('b1', 1) 
        self.n0 = MyDataObject('n0', 2) 
        self.n1 = MyDataObject('n1', 3)
        self.collect = GbmDetectorCollection.from_list([self.b0, self.b1, 
                                                        self.n0, self.n1])
    
    def test_items(self):
        self.assertListEqual(self.collect.items, ['item1', 'item2', 'item3', 'item4'])
    
    def test_types(self):
        self.assertEqual(self.collect.types.__name__, 'MyDataObject')
    
    def test_retrieve_property(self):
        vals = self.collect.value()
        self.assertListEqual(vals, [0, 1, 2, 3])
    
    def test_call_return_function(self):
        val_strs = self.collect.value_as_str()
        self.assertListEqual(val_strs, ['0', '1', '2', '3'])
        
    def test_call_function_all_args(self):
        self.collect.set_value(4)
        self.assertEqual(self.b0.value, 4)
        self.assertEqual(self.b1.value, 4)
        self.assertEqual(self.n0.value, 4)
        self.assertEqual(self.n1.value, 4)

    def test_call_function_bgo_nai_args(self):
        self.collect.set_value(bgo_args=(5,), nai_args=(6,))
        self.assertEqual(self.b0.value, 5)
        self.assertEqual(self.b1.value, 5)
        self.assertEqual(self.n0.value, 6)
        self.assertEqual(self.n1.value, 6)

    def test_call_function_all_kwargs(self):
        res = self.collect.add_int(val=2)
        self.assertEqual(res[0], self.b0.value+2)
        self.assertEqual(res[1], self.b1.value+2)
        self.assertEqual(res[2], self.n0.value+2)
        self.assertEqual(res[3], self.n1.value+2)

    def test_call_function_bgo_nai_kwargs(self):
        res = self.collect.add_int(bgo_kwargs={'val': 3}, nai_kwargs={'val': 4})
        self.assertEqual(res[0], self.b0.value+3)
        self.assertEqual(res[1], self.b1.value+3)
        self.assertEqual(res[2], self.n0.value+4)
        self.assertEqual(res[3], self.n1.value+4)

    def test_remove_and_include(self):
        self.collect.remove('item4')
        self.assertListEqual(self.collect.items, ['item1', 'item2', 'item3'])

        self.collect.include(self.n1, 'n1')
        self.assertListEqual(self.collect.items, ['item1', 'item2', 'item3', 
                                                  'item4'])
    
    def test_no_dets(self):
        collect = GbmDetectorCollection.from_list([0, 1, 2, 3], 
                                                  dets=['b0', 'b1', 'b2', 'b3'])
        self.assertListEqual(collect.items, ['item1', 'item2', 'item3', 'item4'])
        
        collect = GbmDetectorCollection.from_list([0, 1, 2, 3], 
                                                  dets=['b0', 'b1', 'b2', 'b3'],
                                                  names=['o1', 'o2', 'o3', 'o4'])
        self.assertListEqual(collect.items, ['o1', 'o2', 'o3', 'o4'])

    def test_errors(self):
        with self.assertRaises(ValueError):
            collect = GbmDetectorCollection.from_list([0, 1, 2, 3], 
                                                    dets=['b0', 'b1', 'b2'])

        with self.assertRaises(AttributeError):
            collect = GbmDetectorCollection.from_list([0, 1, 2, 3])

        with self.assertRaises(ValueError):
            collect = GbmDetectorCollection.from_list([0, 1, 2, 3], 
                                                    dets=['b0', 'b1', 'b2', 'b3'],
                                                    names=['o1', 'o2', 'o3'])
