
from __future__ import annotations
from nptyping import NDArray, Float64
from typing import Union, TypeVar, NamedTuple, Sequence, Deque, Optional, List, Type, Dict, Any
from types import ModuleType
import xml.etree.ElementTree as ET
from xml.dom.minidom import Document
from xml.dom import minidom
import importlib

from seekr2.libraries.serializer.datatype import data_type

#Generic class type alias
CLASS = TypeVar('CLASS')

#Generic class inst type alias
INSTANCE = TypeVar('INSTANCE')

class SerializerUtils(object):

    def get_object_name(
            self,
            object_: Union[INSTANCE, NamedTuple], 
            ) -> str:
        return object_.__class__.__name__ 

    def get_module_name(
            self,
            class_: CLASS
            ):
        return class_.__module__

    def get_instance_dict(
            self, 
            inst: INSTANCE,
            keys: Optional[bool] = False,
            values: Optional[bool] = False,
            ) -> Union[Dict, List]:
        if keys:
            return list(inst.__dict__)
        elif values:
            return list(inst.__dict__.values())
        return inst.__dict__ 

    def get_namedtuple_dict(
            self, 
            namedtuple_: NamedTuple,
            ):
        return namedtuple_._asdict() 
    
    def get_module(
            self, 
            xml_node: ET.Element,
            ) -> ModuleType:
        if 'module' in xml_node.attrib:
            module_name = self.get_xml_node_attr(xml_node, 'module')
            return importlib.import_module(module_name)
        return importlib.import_module('__main__')
    
    def get_xml_node_attr(
            self, 
            xml_node: ET.Element,
            xml_node_attr_name: str,
            ) -> str:
        return xml_node.attrib[xml_node_attr_name]

    def set_xml_node_attr(
            self, 
            xml_node: ET.Element,
            xml_node_attr_name: str,
            xml_node_attr_value: str,
            ) -> None:
        xml_node.attrib[xml_node_attr_name] = xml_node_attr_value

    def null_adt(
            self,
            size: Optional[int] = None,
            keys: Optional[Union[int, str]] = None,
            key_attr: Optional[str] = None, 
            ) -> Union[List, Dict]:
        if size is not None:
            return self.null_list(size)
        return self.null_dict(keys, key_attr)

    def null_list(
            self,
            size: Optional[int] = None,
            ) -> List:
        return [None for element in range(size)]

    def null_dict(
            self,
            keys: Optional[Any] = None,
            key_attr: Optional[str] = None, 
            ) -> Dict:
        if key_attr is not None:
            keys = [getattr(inst, key_attr) for inst in keys]
            return {key : None for key in keys}
        return {key : None for key in keys}
        

    def class_instance(
            self, 
            xml_node: ET.Element,
            ) -> INSTANCE:
        module = self.get_module(xml_node)
        class_name = None 
        try:
            class_name = self.get_xml_node_attr(xml_node, 'class_name')
        except:
            class_name = self.get_xml_node_attr(xml_node, 'class')
        class_inst = getattr(module, class_name).__new__(
            getattr(module, class_name)
        )
        return class_inst



















