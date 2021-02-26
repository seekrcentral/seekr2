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

from openmmvt.libraries.serializer.datatype import data_type
from openmmvt.libraries.serializer.instanceattrs import InstanceAttrs
from openmmvt.libraries.serializer.strcast import strcast
from openmmvt.libraries.serializer.serializerutils import SerializerUtils
from openmmvt.libraries.serializer.userinputnode import UserInputNode
from openmmvt.libraries.serializer.serialize import Serialize
from openmmvt.libraries.serializer.deserializeuserinput import DeserializeUserInput

#Generic class type alias
CLASS = TypeVar('CLASS')

#Generic class inst type alias
INSTANCE = TypeVar('INSTANCE')

class Deserialize(DeserializeUserInput):

    def set_path(
            self, 
            path: List, 
            index: int,
        ) -> List:
        child_path = copy.deepcopy(path)
        child_path.append(index)
        return child_path

    def deserialize_nested_attr_(
            self, 
            attr, 
            val, 
            path,
            ):
        if len(path) > 1:
            pos = path.popleft()
            return self.deserialize_nested_attr_(
                attr[pos], 
                val, 
                path
            )
        attr[path[0]] = val
        return attr[path[0]]

    def tag_namedtuple(
            self, 
            node,
            ):
        node.tag = (
            node.tag 
            + '___' 
            + self.get_xml_node_attr(node, 'namedtuple_name')
        )
        self.set_xml_node_attr(node, 'type', 'dict')
        return node

    def tag_child_namedtuples(
            self, 
            node,
            ):
        for child in node:
            if self.is_namedtuple(child):
                self.tag_namedtuple(child)

    def is_namedtuple(
            self, 
            node,
            ):
        if self.attr_type(node) == 'namedtuple':
            return True

    def deserialize_nested_attr(
            self, 
            name, 
            val, 
            path,
            ):
        if isinstance(path, list):
            path = deque(path)
        return self.deserialize_nested_attr_(
            self.get_instance_dict(self)[name], 
            val, 
            path
        )

    def deserialize_nested_attr_(
            self, 
            attr, 
            val, 
            path,
            ):
        if len(path) > 1:
            pos = path.popleft()
            return self.deserialize_nested_attr_(
                attr[pos], 
                val, 
                path
            )
        attr[path[0]] = val
        return attr[path[0]]

    def deserialize_list(
            self, 
            node, 
            child = None, 
            path = [],
            ):
        if len(path) != 0:
            null_list = [None for gchild in child]
            name = node.tag
            self.deserialize_nested_attr(name, null_list, path)
            for index, grandchild in enumerate(child):
                grandchild_path = self.set_path(path, index) 
                self.deserialize_child(node, grandchild, grandchild_path)
        else:
            null_list = self.null_adt(size=len(node)) 
            name = node.tag
            setattr(self, name, null_list)
            for index, child in enumerate(node):
                child_path = self.set_path(path, index) 
                self.deserialize_child(node, child, child_path)

    def deserialize_dict(
            self, 
            node, 
            child = None, 
            path = [],
            ):
        if len(path) != 0:
            name = node.tag
            null_dict = self.null_adt(
                keys=child, 
                key_attr='tag'
            )  
            self.deserialize_nested_attr(name, null_dict, path)
            for grandchild in child:
                grandchild_path = self.set_path(path, grandchild.tag) 
                self.deserialize_child(node, grandchild, grandchild_path)
        else:
            null_dict = self.null_adt(
                keys=node, 
                key_attr='tag'
            )  
            name = node.tag
            setattr(self, name, null_dict)
            for child in node:
                child_path = self.set_path(path, child.tag) 
                self.deserialize_child(node, child, child_path)

    def deserialize_namedtuple(
            self, 
            node, 
            child = None, 
            path = [],
            ):
        self.namedtuple_to_dict(node, child, path)
        self.dict_to_namedtuple(self)

    def is_namedtuple_dict(
            self, 
            name,
            ):
        if '___' in name:
            return True

    def namedtuple_name(
            self, 
            dict_name,
            ):
        name = dict_name.split('___')[0]
        return name

    def namedtuple_attrs(
            self, 
            cls,
            ):
        namedtuple_attr_names = []
        non_namedtuple_attr_names = []
        for attr_name in self.get_instance_dict(cls):
            if self.is_namedtuple_dict(attr_name):
                namedtuple_attr_names.append(attr_name)
            else:
                non_namedtuple_attr_names.append(attr_name)
        return namedtuple_attr_names, non_namedtuple_attr_names
            

    def namedtuple_to_dict(
            self, 
            node, 
            child = None, 
            path = [],
            ):
        name = node.tag
        if self.is_namedtuple(node):
            name = self.tag_namedtuple(node) 
        if len(path) != 0:
            self.tag_child_namedtuples(child)
            null_dict = self.null_adt(
                keys=child, 
                key_attr='tag'
            )  
            self.deserialize_nested_attr(name, null_dict, path)
            for grandchild in child:
                grandchild_path = self.set_path(path, grandchild.tag) 
                self.deserialize_child(node, grandchild, grandchild_path)
        else:
            name = node.tag
            self.tag_child_namedtuples(node)
            null_dict = self.null_adt(
                keys=node, 
                key_attr='tag'
            )  
            setattr(self, name, null_dict)
            for child in node:
                child_path = self.set_path(path, child.tag) 
                self.deserialize_child(node, child, child_path)

    def dict_to_namedtuple(
            self, 
            mapping, 
            name = '',
            ):
        if data_type(mapping) == 'instance':
            namedtuple_attrs, not_namedtuple_attrs \
                = self.namedtuple_attrs(mapping)  
            for key in namedtuple_attrs:
                setattr(
                    mapping, 
                    key.split('___')[0], 
                    self.dict_to_namedtuple(
                        getattr(mapping, key), 
                        key
                    )
                )
        if isinstance(mapping, dict):
            for key in list(mapping):
                if '___' in key:
                    mapping[key.split('___')[0]] \
                        = self.dict_to_namedtuple(mapping[key], key)
                else:
                    mapping[key] \
                        = self.dict_to_namedtuple(mapping[key], key)
            return self.namedtuple_from_mapping(mapping, name)
        return mapping

    def namedtuple_from_mapping(
            self, 
            mapping, 
            name,
            ):
        if self.is_namedtuple_dict(name):
            mapping_keys = [key.split('___')[0] for key in list(mapping)]
            mapping_values = []
            for key in mapping_keys:
                mapping_values.append(mapping[key])
            mapping_ = dict(zip(mapping_keys, mapping_values))
            this_namedtuple_maker = namedtuple(name.split('___')[1], mapping_)
            return this_namedtuple_maker(**mapping_)
        else:
            return mapping
         

    def deserialize_instance(
            self, 
            node, 
            child = None, 
            path = [],
            ):
        if len(path) != 0:
            name = node.tag 
            class_inst = self.class_instance(child)
            self.deserialize_nested_attr(name, class_inst, path)
            class_inst._deserialize(child)
        else:
            class_inst = self.class_instance(node)
            setattr(self, node.tag, class_inst)
            class_inst._deserialize(node)
        
    def attr_type(
            self, 
            node,
            ):
        if not hasattr(node, 'attrib'):
            return 'NoneType'
        return node.attrib['type']

    def deserialize_child(
            self, 
            node, 
            child = None, 
            path = [],
            ):
        deserialize_attr_dict = {
            'int' : self.deserialize_int,
            'float' : self.deserialize_float,
            'str' : self.deserialize_str,
            'bool' : self.deserialize_bool,
            'instance' : self.deserialize_instance,
            #Added for legacy support
            'class' : self.deserialize_instance, 
            'list' : self.deserialize_list,
            'dict' : self.deserialize_dict,
            'namedtuple' : self.deserialize_namedtuple,
            'NoneType' : self.deserialize_nonetype
        }
        if len(path) != 0:
            typ = self.attr_type(child)
            return deserialize_attr_dict[typ](node, child, path)
        typ = self.attr_type(node)
        return deserialize_attr_dict[typ](node, child, path)

    def deserialize_int(
            self, 
            node, 
            child = None, 
            path = [],
            ):
        if path:
            name = node.tag
            val = int(child.text)
            return self.deserialize_nested_attr(name, val, path)
        name = node.tag
        val = int(node.text)
        setattr(self, name, val)

    def deserialize_nonetype(
            self, 
            node, 
            child = None, 
            path = [],
            ):
        if path:
            name = node.tag
            val = None 
            return self.deserialize_nested_attr(name , val, path)
        name = node.tag
        val = None 
        setattr(self, name, val)

    def deserialize_float(
            self, 
            node, 
            child = None, 
            path = [],
            ):
        if path:
            name = node.tag
            val = float(child.text)
            return self.deserialize_nested_attr(name , val, path)
        name = node.tag
        val = float(node.text)
        setattr(self, name, val)

    def deserialize_bool(
            self, 
            node, 
            child = None, 
            path = [],
            ):
        if path:
            name = node.tag
            val = strcast(child.text)
            return self.deserialize_nested_attr(name, val, path)
        name = node.tag
        val = strcast(node.text)
        setattr(self, name, val)

    def set_str(
            self, 
            s,
            ):
        if s == 'None':
            return ''
        return s
        
    def deserialize_str(
            self, 
            node, 
            child = None, 
            path = [],
            ):
        if path:
            name = node.tag
            val = self.set_str(str(child.text))
            return self.deserialize_nested_attr(name, val, path)
        name = node.tag
        val = self.set_str(str(node.text))
        setattr(self, name, val)

    def _deserialize(
            self, 
            root, 
            path = [],
            ):
        for child in root:
            if 'type' in child.attrib:
                self.deserialize_child(child)
            else:
                setattr(self, child.tag, None)
                            
    def deserialize(
            self, 
            root, 
            user_input = False,
            ):
        xml_root = ET.parse(root).getroot()
        if user_input:
            inst_attrs = InstanceAttrs(self)
            empty_str_attrs = []
            for attr in inst_attrs.default_instance_attrs():
                if inst_attrs.is_empty_str(attr):
                    empty_str_attrs.append(attr)
            self.deserialize_user_input(
                root, 
                xml_root, 
                empty_str_attrs, 
                inst_attrs
            )
            return
        self._deserialize(xml_root)
                    

'''
-----------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------
'''

class Serializer(Deserialize, Serialize):
    pass

