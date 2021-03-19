from __future__ import annotations
from nptyping import NDArray, Float64
from typing import Union, TypeVar, NamedTuple, Sequence, Deque, Optional, List, Type, Dict, Any
from types import ModuleType
import copy
import inspect
import ast
import textwrap
import importlib
from seekr2.libraries.serializer.serializerutils import SerializerUtils
from functools import lru_cache as lru


#Generic class type alias
CLASS = TypeVar('CLASS')

#Generic class inst type alias
INSTANCE = TypeVar('INSTANCE')

class InstanceAttrs(SerializerUtils):

    def __init__(
            self,
            inst: INSTANCE,
            ) -> None:
        self.inst = inst 

    def init_src_code(self) -> str:
        return textwrap.dedent(
            inspect.getsource(self.inst.__init__))

    def instance_attrs_ast(
            self,
            init_src_code: str,
            ) -> List[ast.Assign]:
        inst_attrs_ast = ast.parse(init_src_code).body[0].body
        return [
            node for node in inst_attrs_ast 
            if isinstance(node, ast.Assign)
        ]


    def default_instance_attr_names(
            self,
            inst_attrs_ast: List[ast.Assign],
            ) -> List[str]:
        return [
            node.targets[0].attr for node in inst_attrs_ast
            if isinstance(node, ast.Assign)
        ]

    def default_instance_attr_values(
            self,
            inst_attrs_ast: List[ast.Assign],
            ):
        for node in inst_attrs_ast:
            value = self.get_instance_dict(node.value, values=True)[0]
            yield value

    def default_instance_attrs(self) -> Dict[str, Any]:
        init_src_code = self.init_src_code() 
        inst_attrs_ast = self.instance_attrs_ast(init_src_code) 
        default_inst_attr_names = \
            self.default_instance_attr_names(inst_attrs_ast)
        default_inst_attr_values = \
            list(self.default_instance_attr_values(inst_attrs_ast))
        return dict(zip(
            default_inst_attr_names, default_inst_attr_values
            ))

    def is_empty_str(
            self,
            default_inst_attr: str,
            ) -> Union[bool, None]:
        default_inst_attrs = self.default_instance_attrs()
        if default_inst_attrs[default_inst_attr] == '':
            return True

    def init_args(self):
        class_ = type(self.inst)
        signature = inspect.signature(self.inst.__init__)
        return [
            arg_name for arg_name in list(signature.parameters)
        ]

    def instance_attrs(self):
        return copy.deepcopy(self.get_instance_dict(self.inst))

    def map_instance_attrs_to_init_args(self):
        inst_attr_names = list(self.instance_attrs())
        args = []
    
        for arg in self.init_args():
            if arg in inst_attr_names:
                args.append(getattr(self.inst, arg))
            else:
                return
        return tuple(args)

    def call_init__(
            self,
            args, 
            ):
        inst_attrs = self.instance_attrs()
        self.inst.__init__(*args)
        for attr_name, attr_value in inst_attrs.items():
            setattr(self.inst, attr_name, attr_value)  
        return

    def call_init(self):
        args = self.map_instance_attrs_to_init_args()
        if isinstance(args, tuple):
            self.call_init__(args)

