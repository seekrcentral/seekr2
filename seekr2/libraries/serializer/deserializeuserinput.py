from __future__ import annotations
from nptyping import NDArray, Float64
from typing import Union, TypeVar, NamedTuple, Sequence, Deque, Optional, List, Type, Dict, Any
from types import ModuleType
from collections import namedtuple
import copy
import inspect
import xml.etree.ElementTree as ET
from xml.dom.minidom import Document
from xml.dom import minidom
import importlib
from collections import deque

from seekr2.libraries.serializer.datatype import data_type
from seekr2.libraries.serializer.instanceattrs import InstanceAttrs
from seekr2.libraries.serializer.strcast import strcast
from seekr2.libraries.serializer.serializerutils import SerializerUtils
from seekr2.libraries.serializer.userinputnode import UserInputNode
from seekr2.libraries.serializer.serialize import Serialize


#Generic class type alias
CLASS = TypeVar('CLASS')

#Generic class inst type alias
INSTANCE = TypeVar('INSTANCE')

#Add documentation...otherwise good



class DeserializeUserInput(SerializerUtils):

    def init_args(
            self, 
            inst: INSTANCE,
        ) -> None:
        class_ = type(inst)
        signature = inspect.signature(class_)
        init_arg_names = [
            arg_name for arg_name in list(signature.parameters)
        ]
        inst_attr_names = self.get_inst_dict(keys=True) 
        

    def instance_attr(
            self, 
            node: UserInputNode, 
            xml_node: ET.Element, 
            inst_attrs: InstanceAttrs,
        ) -> None:
        module = importlib.import_module('__main__')
        inst = self.class_instance(xml_node)
        new_inst_attrs = InstanceAttrs(inst)
        empty_str_attrs = []
        for attr in new_inst_attrs.default_instance_attrs():
            if new_inst_attrs.is_empty_str(attr):
                empty_str_attrs.append(attr)
        inst_attrs.call_init()
        setattr(self, node.name, inst)
        inst.__deserialize_user_input(
            node, 
            xml_node, 
            empty_str_attrs, 
            new_inst_attrs
        )

    def instance_list_attr(
            self, 
            node: UserInputNode, 
            xml_node: ET.Element, 
            inst_attrs: InstanceAttrs,
        ) -> None:
        num_children = len(xml_node.getchildren())
        null_list = self.null_adt(num_children)
        inst_list = []
        children = []
        inst_list_empty_str_attrs = []
        inst_list_attrs = []
        for index in range(len(null_list)):
            inst = self.class_instance(xml_node[index])
            new_inst_attrs = InstanceAttrs(inst)
            empty_str_attrs = []
            for attr in new_inst_attrs.default_instance_attrs():
                if new_inst_attrs.is_empty_str(attr):
                    empty_str_attrs.append(attr)
            inst_list.append(inst)
            inst_list_attrs.append(new_inst_attrs)
            inst_list_empty_str_attrs.append(empty_str_attrs)
            children.append(node.add_child(xml_node[index]))
        inst_attrs.call_init()
        setattr(self, node.name, inst_list) 
        for index in range(len(children)):
            inst_list[index].__deserialize_user_input(
                children[index], 
                xml_node[index], 
                inst_list_empty_str_attrs[index], 
                inst_list_attrs[index]
            )

    def non_instance_attr(
            self, 
            node: UserInputNode, 
            xml_node: ET.Element, 
            empty_str_attrs: List, 
            inst_attrs: InstanceAttrs,
        ) -> None:
        if (type(xml_node.text) == type(None)
                and xml_node.tag in empty_str_attrs):
            setattr(self, node.name, '')
            inst_attrs.call_init()
        else:
            setattr(
                self, 
                node.name, 
                strcast(xml_node.text)
            )
            inst_attrs.call_init()

    def __deserialize_user_input(
            self, 
            node: UserInputNode, 
            xml_node: ET.Element, 
            empty_str_attrs: List, 
            inst_attrs: InstanceAttrs,
            ) -> UserInputNode:
        for xml_child in node.xml_node:
            child = node.add_child(xml_child)
            #Added "or child.type == 'class'" for legacy support
            if child.type_ == 'instance' or child.type_ == 'class': 
                self.instance_attr(child, xml_child, inst_attrs)
            elif child.type_ == 'instance_list':
                self.instance_list_attr(child, xml_child, inst_attrs)
            elif child.type_ == 'non_instance':
                self.non_instance_attr(
                    child, 
                    xml_child, 
                    empty_str_attrs, 
                    inst_attrs,
                )
        return node
            

    def deserialize_user_input(
            self, 
            node: UserInputNode, 
            xml_root: ET.Element, 
            empty_str_attrs: List, 
            inst_attrs: InstanceAttrs,
        ) -> str:
        name = xml_root.tag
        root = UserInputNode(xml_root, name, 'instance')
        self.__deserialize_user_input(
            root, 
            xml_root, 
            empty_str_attrs, 
            inst_attrs, 
        )
        return root
'''
-----------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------
'''
