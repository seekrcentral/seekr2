from __future__ import annotations
from typing import Union, TypeVar, NamedTuple, Sequence, Deque, Optional, List, Type, Dict, Any
from types import ModuleType

from nptyping import NDArray, Float64

from openmmvt.libraries.serializer.strcast import strcast



CLASS = TypeVar('CLASS')

#Generic class instance type alias
INSTANCE = TypeVar('INSTANCE')

class UserInputNode(object):
    
    def __init__(self,
            xml_node: ET.Element, 
            name: str, 
            type_: str, 
            parent: Optional[UserInputNode] = None,
            ) -> None:
        self.xml_node = xml_node
        self.name = name
        self.type_ = type_ 
        self.parent = parent
        self.val = None

    def data_type(
            self, 
            xml_node,
            ):
        if not xml_node.getchildren():
            return 'non_instance'
        else:
            class_children = []
            for child in xml_node.getchildren():
                if getattr(child, 'attrib'):
                    class_children.append(child)
            if len(class_children) == len(xml_node.getchildren()):
                return 'instance_list'
            else:
                return 'instance'

    def get_name(
            self, 
            xml_node,
            ):
        return xml_node.tag 

    def add_child(
            self, 
            xml_node,
            ):
        name = self.get_name(xml_node)
        type_ = self.data_type(xml_node)
        child = UserInputNode(
            xml_node, 
            name, type_, 
            self
        )
        if type_ == 'pdt':
            child.data = strcast(xml_node.text)
        setattr(self, name, child)
        return getattr(self, name)
