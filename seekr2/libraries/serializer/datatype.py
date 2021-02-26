from __future__ import annotations
from collections.abc import Iterable
import re
from nptyping import NDArray, Float64
from typing import Union, TypeVar, NamedTuple, Sequence, Deque, Optional, List, Type, Dict, Any
from types import ModuleType

class CheckType(object):
    """
    CheckType is a class that is used to determine if a specific data object 
    is a user defined class instance, namedtuple, list, or dictionary. 
    
    """

    @classmethod
    def is_instance(
            cls,
            data: Any,
            ) -> bool:
        if hasattr(data, '__dict__'):
            if hasattr(data, '__module__'):
                return True

    @classmethod
    def is_namedtuple(
            cls,
            data: Any,
            ) -> bool:
        typ = type(data)
        bases = typ.__bases__
        if len(bases) != 1 or bases[0] != tuple:
            return False
        fields = getattr(typ, '_fields', None)
        if not isinstance(fields, tuple):
            return False
        return all(type(ele) == str for ele in fields)

    @classmethod
    def is_dict(
            cls,
            data: Any,
            ) -> bool:
        if isinstance(data, dict):
            return True 

    @classmethod
    def is_list(
            cls,
            data: Any,
            ) -> bool:
        if isinstance(data, list):
            return True 

def data_type(
        data: Any,
        ) -> str: 
    if CheckType.is_instance(data):
        return 'instance'
    elif CheckType.is_namedtuple(data):
        return 'namedtuple'
    elif CheckType.is_dict(data):
        return 'dict'
    elif CheckType.is_list(data):
        return 'list'


