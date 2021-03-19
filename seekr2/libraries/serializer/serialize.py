
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

#Generic class type alias
CLASS = TypeVar('CLASS')

#Generic class inst type alias
INSTANCE = TypeVar('INSTANCE')


class Serialize(SerializerUtils):

    def _delete_class_vars(self) -> None:
        if hasattr(self, 'doc'):
            del self.__class__.doc
            del self.__class__.rootName
            del self.__class__.root

    def _set_class_vars(self) -> None:
        self.__class__.doc = Document()
        self.__class__.rootName = self.__class__.__name__
        self.__class__.root = self.doc.createElement(self.rootName)
        self.doc.appendChild(self.root)

    def data_type(
            self, 
            data: Any,
            ) -> str:
        return data.__class__.__name__

    def is_abstract_data_type(
            self, 
            data: Any,
            ) -> Union[None, bool]:
        data_type = self.data_type(data)
        pdts = set([
            'int', 'float', 'bool', 
            'str', 'float64', 'NoneType'
        ])
        if data_type in pdts:
            return
        return True

    def pdt_str_cast(
            self, 
            data: Union[
                int, str, float,
                bool, Float64, None,
            ]
            ) -> str:
        return str(data)

    def set_xml_attribs(
            self, 
            parent: ET.Element, 
            **kwargs,
            ) -> None:
        for attrib_name, attrib_value in kwargs.items():
            parent.setAttribute(attrib_name, attrib_value)
    
    def set_list_tag(
            self, 
            index: int, 
            minidom_doc: Document, 
            parent: ET.Element,
            ) -> ET.Element:
        return minidom_doc.createElement((
            parent.tagName
            + '_e'
            + str(index)
        ))

    def set_dict_tag(
            self, 
            key, 
            minidom_doc, 
            parent
        ) -> ET.Element:
        tag = None
        if isinstance(key, int):
            tag = minidom_doc.createElement((
                'key_'
                + str(key)
            ))
        else:
            tag = self.doc.createElement(key)
        return tag
            
    def serialize_pdt(
            self, 
            data, 
            minidom_doc, 
            parent,
        ) -> None: 
        xml_str = self.pdt_str_cast(data)
        tag = minidom_doc.createTextNode(xml_str)
        parent.setAttribute('type', self.data_type(data))
        parent.appendChild(tag)

    def serialize_struct(
            self, 
            struct, 
            minidom_doc, 
            parent,
        ) -> None:
        typ = data_type(struct)
        if typ == 'instance':
            self.serialize_inst(struct, minidom_doc, parent)
        elif typ == 'namedtuple':
            self.serialize_namedtuple(struct, minidom_doc, parent)
        elif typ == 'list':
            self.serialize_list(struct, minidom_doc, parent)
        elif typ == 'dict':
            self.serialize_dict(struct, minidom_doc, parent)

    def serialize_inst(
            self, 
            cls, 
            minidom_doc, 
            parent,
        ) -> None:
        for attr_name, attr_value in self.get_instance_dict(cls).items():
            attr_value = self.get_instance_dict(cls)[attr_name]
            tag = minidom_doc.createElement(attr_name)
            self.set_xml_attribs(
                parent, 
                type = 'instance', 
                class_name = self.get_object_name(cls),
                module = self.get_module_name(cls)
            ) 
            parent.appendChild(tag)
            self._serialize(attr_value, tag)

    def serialize_namedtuple(
            self, 
            namedtuple_: NamedTuple, 
            minidom_doc: Document, 
            parent: ET.Element,
        ) -> None:
        for field_key, field_val in namedtuple_._asdict().items():
            if field_key.startswith('__'):
                continue
            tag = minidom_doc.createElement(field_key)
            self.set_xml_attribs(
                parent, 
                type = 'namedtuple',
                namedtuple_name = self.get_object_name(namedtuple_) 
            )
            parent.appendChild(tag)
            self._serialize(field_val, tag)

    def serialize_list(
            self, 
            list_: List, 
            minidom_doc: Document, 
            parent: ET.Element,
        ) -> None:
        for i, e in enumerate(list_):
            tag = self.set_list_tag(i, minidom_doc, parent) 
            parent.appendChild(tag)
            self.set_xml_attribs(parent, type = 'list')
            parent.setAttribute('type', 'list')
            self._serialize(e, tag)

    def serialize_dict(
            self, 
            dict_: Dict, 
            minidom_doc: Document, 
            parent: ET.Element,
        ) -> None:
        for key in dict_:
            tag = self.set_dict_tag(key, minidom_doc, parent)
            parent.appendChild(tag)
            self.set_xml_attribs(parent, type = 'dict')
            self._serialize(dict_[key], tag)

    def add_child(
            self, 
            data: Any, 
            minidom_doc: Document, 
            parent: ET.Element,
        ) -> None:
        if self.data_type(data) == 'NoneType':
            return
        elif self.is_abstract_data_type(data):
            self.serialize_struct(data, minidom_doc, parent)
        else:
            self.serialize_pdt(data, minidom_doc, parent)
        

    def _serialize(
            self, 
            value: Any, 
            parent: ET.Element,
        ) -> None:
        self.add_child(value, self.doc, parent)

    def serialize(
            self, 
            xml_filename: Optional[str] = None,
        ) -> None:
        self._set_class_vars()
        parent = self.root
        self._serialize(self, parent)
        if xml_filename is not None:
            self._write_xml(xml_filename)
        self._delete_class_vars()

    def _write_xml(
            self, 
            xml_filename: str,
        ) -> None:    
        xmlstr = self.doc.toprettyxml(indent="    ")
        xml_file = open(xml_filename, 'w')
        xml_file.write(xmlstr)
        xml_file.close()

