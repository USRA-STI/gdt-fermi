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
from gdt.core.collection import DataCollection

__all__ = ['GbmDetectorCollection']

class GbmDetectorCollection(DataCollection):
    """A container for a collection of GBM-specific data objects, such as a 
    collection of :class:`~.phaii.GbmPhaii` objects from different detectors.
    
    The special behavior of this class is to provide a way to interact with
    a collection of detector data that may contain a mix of different *types*
    of detectors.  For example, many times we want a collection of GBM NaI
    and GBM BGO detectors.  These detectors have very different energy ranges,
    and so may require different inputs for a variety of functions.  This
    collection allows one to specify the different arguments for NaI and BGO
    data without having to implement many ugly and space-wasting loops and
    ``if...else`` blocks.
    
    In addition to the GbmDetectorCollection methods, all of the individual 
    object attributes and methods are exposed, and they become methods of the 
    DataCollection.  Note that individual object attributes become *methods* 
    i.e. if you have an item attribute called item.name, then the corresponding 
    DataCollection method would be item.name().    
    """
    def __init__(self):
        super().__init__()
        self._dets = []

    @classmethod
    def from_list(cls, data_list, names=None, dets=None):
        """Given a list of objects and optionally a list of corresponding names
        and corresponding detector names, create a new GbmDetectorCollection. 
        
        Args:
            data_list (list of :obj:`objects`): 
                The list of objects to be in the collection
            names (list of :obj:`str`, optional):  
                The list of corresponding names to the objects.  If not set, 
                will try to retrieve a name from object.filename (assuming it's 
                a data object). If that fails, each item will be named 
                ambiguously 'item1', 'item2', etc.
            dets (list of :obj:`str`, optional): 
                The detector names for each object. If not set, will try to
                retrieve from the object.detector attribute.  If that attribute 
                doesn't exist, an error will be raised, and the user will need 
                to specify this list.
        
        Returns                
            :py:class:`GbmDetectorCollection`: The newly created collection
        """
        obj = cls()

        # set the detector names
        if dets is not None:
            if len(dets) != len(data_list):
                raise ValueError(
                    'Detector list must be same size as data list')
        else:
            try:
                dets = [data_item.detector for data_item in data_list]
            except:
                raise AttributeError('Cannot find detector information. '
                                     'Need to manually set')
        # set the names
        if names is not None:
            if len(names) != len(data_list):
                raise ValueError('Names list must be same size as data list')
        else:
            names = [None] * len(data_list)

        # include the objects
        [obj.include(data_item, det, name=name) for (data_item, det, name)
         in zip(data_list, dets, names)]

        return obj

    def remove(self, item_name):
        """Remove an object from the collection given the name 
        
        Args:
            item_name (str): The name of the item to remove
        """
        index = [item == item_name for item in self.items].index(True)
        self._dets.pop(index)
        self._data_dict.pop(item_name)

    def include(self, data_item, det, name=None):
        """Insert an object into the GbmDetectorCollection.  The first item 
        inserted will set the immutable type.
        
        Args:
            data_item (:obj:`object`): A data object to include
            det (str): The corresponding detector for the item
            name (str, optional): 
                An optional corresponding name.  If not set, will try to 
                retrieve a name from object.filename (assuming it's a data 
                object). If that fails, each item will be named ambiguously 
                'item1', 'item2', etc.
        """
        super().include(data_item, name=name)
        self._dets.append(det)

    def _method_call(self, method_name, *args, nai_args=(), nai_kwargs=None,
                     bgo_args=(), bgo_kwargs=None, **kwargs):
        """This is the wrapper for the attribute and method calls.  Applies 
        method_name over all items in the GbmDetectorCollection.

        Args:
            method_name (str): The name of the method or attribute
            *args: Additional arguments to be passed to the method
            nai_args: Arguments to be applied only to the NaI objects
            bgo_args: Arguments to be applied only to the BGO objects
            nai_kwargs: Keywords to be applied only to the NaI objects
            bgo_kwargs: Keywords to be applied only to the BGO objects
            **kwargs: Additional keyword arguments to be passed to the 
                      method. Will be applied to both NaI and BGO objects 
                      and will be appended to any existing keywords from 
                      nai_kwargs or bgo_kwargs
        
        Returns:       
        None or list: If not None, will return the results from all objects 
                      in the list     
        """
        if nai_kwargs is None:
            nai_kwargs = {}
        if bgo_kwargs is None:
            bgo_kwargs = {}

        if len(args) > 0:
            nai_args = args
            bgo_args = args
        if len(kwargs) > 0:
            nai_kwargs.update(kwargs)
            bgo_kwargs.update(kwargs)

        # get the attributes/methods for each item
        refs = [getattr(obj, method_name) for obj in self._data_dict.values()]

        # if method_name is a method, then it will be callable
        if callable(refs[0]):
            res = []
            for obj, det in zip(self._data_dict.values(), self._dets):
                if 'n' in det:
                    our_args = nai_args
                    our_kwargs = nai_kwargs
                elif 'b' in det:
                    our_args = bgo_args
                    our_kwargs = bgo_kwargs
                res.append(getattr(obj, method_name)(*our_args, **our_kwargs))

        # otherwise, method_name will not be callable if it is an attribute
        else:
            # we are setting an attribute    
            if len(nai_args) != 0 or len(bgo_args) != 0:
                res = []
                for obj, det in zip(self._data_dict.values(), self._dets):
                    if 'n' in obj.detector:
                        our_args = nai_args
                    elif 'b' in obj.detector:
                        our_args = bgo_args
                    res.append(setattr(obj, method_name, *args))
            # we are retrieving an attribute
            else:
                res = refs

        if res[0] is not None:
            return res
