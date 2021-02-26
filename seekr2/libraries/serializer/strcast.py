import numpy as np
from collections import deque
import re
import sys


class ListStrCastNode(object): 

    def __init__(
            self,
            start,
            end,
            parent = None,
            path = [],
            ):
        self.start = start
        self.end = end
        self.parent = parent
        self.utils = StrCastUtils()
        self.path = path 

    def top_child(
            self,
            children,
            ):
        return list(children)[-1] 

    def pop_child(
            self,
            children,
            pos,
            ):
        return children.pop(pos) 

    def is_child(
            self,
            start1,
            end1,
            start2,
            end2,
            ):
        if start1 > start2 and end1 < end2:
            return True
        return

    def is_leaf(
            self,
            children,
            pos,
            ):
        if self.top_child(children) == pos:
            return True
        return

    def is_chars_leaf(
            self,
            chars_dict,
            pos,
            ):
        if self.top_child(chars_dict) == pos:
            return True
        return

    def get_children(
            self,
            internal_nodes = True,
            leaves = True,
            ):
        if internal_nodes and leaves:
            return [
                child for name, child in self.__dict__.items() 
                    if 'child' in name or 'leaf' in name
            ]
        elif internal_nodes and not leaves:
            return [
                child for name, child in self.__dict__.items() 
                    if 'child' in name
            ]
        elif leaves and not internal_nodes:
            return [
                child for name, child in self.__dict__.items() 
                    if 'leaf' in name
            ]
        return

    def child_names(
            self, 
            node,
            ):
        return [
            name for name in self.__dict__ 
                if 'child' in name
        ]

    def leaf_names(
            self,
            node,
            ):
        return [
            name for name in self.__dict__
                if 'leaf' in name
        ]

    def child_path(self):
        path = [] 
        if self.path:
            path = [pos for pos in self.path] 
        path.append(len(self.get_children(self)))
        return path

    def add_child(
            self, 
            start,
            end,
            ):
        if self.child_names(self):
            name = (
                'child'
                + str(int(re.findall(
                    r'\d+', 
                    self.child_names(self)[-1])[0]) 
                + 1)
            )
        else:
            name = 'child0'
        setattr(
            self, name, ListStrCastNode(
                start, end, self,
                self.child_path()
            )
        ) 
        return getattr(self, name) 

    def add_leaf(
            self,
            start,
            end, 
            s,
            chars = False
            ):
        name = None
        if self.leaf_names(self):
            name = (
                'leaf'
                + str(int(re.findall(
                    r'\d+', 
                    self.leaf_names(self)[-1])[0])
                + 1)
            )
        else:
            name = 'leaf0'
        leaf = self.utils.cleanup_str(
            s, start + 1,
            end, 'list'
        )
        if chars:
            setattr(self, name, tuple(leaf))
        else:
            setattr(self, name, leaf)
        return getattr(self, name) 


class StrCastUtils(object):

    def xml_str_cast(
            self,
            val,
            ):
        if type(val.text) == type(None):
            return   
        else:
            if val.text.lower() == 'true':
                return True 
            elif val.text.lower() == 'false':
                return False 
            elif '[' in val.text:
                cast = StrCast()
                return cast.str_to_list(val)
            elif val.text == '.' or re.search('[a-zA-Z]', val.text):
                return self.str_cast_(val.text, 'str')
            elif re.search('[\.]', val.text):
                return self.str_cast_(val.text, 'float')
            else:
                return self.str_cast_(val.text, 'int')

    def std_str_cast(
            self,
            val,
            ):
        if type(val) == type(None):
            return   
        else:
            if val.lower() == 'true':
                return True 
            elif val.lower() == 'false':
                return False 
            elif '[' in val:
                cast = StrCast()
                return cast.str_to_list(val)
            elif val == '.' or re.search('[a-zA-Z]', val):
                return self.__str_cast(val, 'str')
            elif re.search('[\.]', val):
                return self.__str_cast(val, 'float')
            else:
                return self.__str_cast(val, 'int')

    def str_cast(
            self, 
            val,
            ):
        if hasattr(val, 'text'):
            return self.xml_str_cast(val)
        else:
            return self.std_str_cast(val)

    def __str_cast(
            self, 
            str_val,
            type_str,
            ):
        type_lookup = {
            'int' : int,
            'float' : float,
            'str' : str
        } 
        return type_lookup[type_str](str_val)

    def except_the_kitchen_sink_regex(
            self, 
            str_struct
            ):
        struct = re.findall(
            r'\S*[A-Za-z]\S*[0-9]\S*|'
            '[0-9]?\S*[A-Za-z]\S*|'
            '\.?\d+\.?\d*|\d+',
            str_struct
        )
        struct = [char.replace(',', '') for char in struct]
        return struct

    def str_cast_list_elements(
            self,
            str_list,
            ):
        return [
            self.str_cast(char) for char 
            in self.except_the_kitchen_sink_regex(str_list)
        ]

    def str_cast_dict_elements(
            self,
            str_dict,
            ): #Come back to this later
        return  

#get rid of *args
    def str_cast_struct_elements(
            self,
            str_struct,
            *args,
            ):
        str_cast_struct_dict = {
            'list' : self.str_cast_list_elements(str_struct), 
            'dict' : self.str_cast_dict_elements(str_struct),
        }
        struct_type = args[0]
        return str_cast_struct_dict[struct_type]

#get rid of *args
    def find_matching(
            self,
            s,
            *args,
            ):
        matching_dict = {}
        stack = []
        start, end = args
        for pos, i in enumerate(s):
            if i == start:
                stack.append(pos)
            elif i == end:
                matching_dict[stack.pop()] = pos
        return matching_dict

    def num_lists(
            self,
            s,
            ):
        return self.num_structs(s, 'list')

    def num_dicts(
            self,
            s,
            ):
        return self.num_structs(s, 'dict')

    def num_structs(
            self,
            s,
            struct_type,
            ):
        num_structs_dict = { 
            'list' : self.find_matching(s, '[', ']'), 
            'dict' : self.find_matching(s, '{', '}')
        }
        return len(num_structs_dict[struct_type])

    def nested_list_insert(
            self,
            lst,
            val,
            path,
            ):
        if isinstance(path, list):
            path = deque(path)
        self.__nested_list_insert(lst, val, path)
        return lst

    def __nested_list_insert(
            self,
            lst,
            val,
            path,
            ):
        if path:
            pos = path.popleft()
            return self.__nested_list_insert(lst[pos], val, path)
        lst.append(val)

    def sort_dict(
            self,
            d,
            *args,
            **kwargs,
            ):
        args = list(args)
        for index, arg in enumerate(args):
            if arg == 'key':
                args[index] = 0
            if arg == 'value':
                args[index] = 1
        kwargs_keys = [key for key, value in kwargs.items()]
        return {
            key : value for key, value in sorted(
                d.items(), key = lambda item : item[args[0]], 
                reverse = kwargs[kwargs_keys[0]]
            )
        }

    def cleanup_str(
            self, 
            s, 
            start, 
            end, 
            *args,
            ):
        return self.str_cast_struct_elements(s[start : end], *args)

    def chars(
            self, 
            s, 
            *args,
            ):
        chars_dict = {}
        prev_bracket = {}
        for index, char in enumerate(s):
            if char == args[0]:
                if prev_bracket:
                    start, end = list(prev_bracket)[-1], index
                    if prev_bracket[start] == args[0]:
                        if (end - start) != 1:
                            chars_dict[start] = self.cleanup_str(
                                s, start, end,
                                args[2]
                            )
                    elif prev_bracket[start] == args[1]:
                        if (end - start) != 3:
                            chars_dict[start] = self.cleanup_str(
                                    s, start, end,
                                    args[2]
                                )
                prev_bracket = {index : args[0]} 
            elif char == args[1]:
                start, end = list(prev_bracket)[-1], index
                if prev_bracket[start] == args[1]:
                    if (end - start) != 1:
                        chars_dict[start] = self.cleanup_str(
                            s, start, end,
                            args[2]
                        )
                prev_bracket = {index : args[1]} 
        return chars_dict

    
class StrCast(object):

    def __init__(self):
        self.utils = StrCastUtils()

    def adjust_recursion_limit(
            self,
            s,
            ):
        if self.utils.num_lists(s) > 100:
            sys.setrecursionlimit(int(10e8))

    def start_end(
            self, 
            q_s, 
            pos = None,
            ):
        if not pos:
            start = list(q_s[-1])[0]
            end = q_s[-1][start]
            return start, end 
        start = list(q_s[pos])[0]
        end = q_s[pos][start]
        return start, end 
        
    def list_init_queue(
            self, 
            s,
            ):
        return deque([
            {start : end} for start, end in self.utils.sort_dict(
                self.utils.find_matching(s, '[', ']'),
                'key', reverse = True).items()
        ])

    def list_init_children(
            self, 
            s,
            ):
        return self.utils.sort_dict(
            self.utils.find_matching(s, '[', ']'),
            'value', reverse = True
        )

    def list_init_stack(
            self,
            s,
            init_queue = None,
            ):
        if init_queue:
            return deque([init_queue.pop()])
        else:
            return deque()

    def list_str_tree(
            self,
            s,
            **qscpn,
            ):
        if not qscpn:
            self.adjust_recursion_limit(s)
            qscpn['queue'] = self.list_init_queue(s)
            if len(qscpn['queue']) == 1:
                qs, qe = self.start_end(qscpn['queue']) 
                node = ListStrCastNode(qs, qe)
                node.add_leaf(qs, qe, s)
                return node
            qscpn['stack'] = self.list_init_stack(s, qscpn['queue'])
            qscpn['children'] = self.list_init_children(s) 
            qscpn['pos'] =  self.start_end(qscpn['queue'])[0]
            qscpn['node'] = ListStrCastNode(*self.start_end(qscpn['stack']))
            qscpn['node'].pop_child(qscpn['children'], 0)
        queue, stack, children, pos, node = qscpn.values()
        qs, qe = self.start_end(queue) 
        ss, se = self.start_end(stack) 
        if len(queue) > 1:
            if se > qe:
                if node.is_leaf(children, pos):
                    node.pop_child(children, qs)
                    node.add_leaf(qs, qe, s)
                    queue.pop()
                    pos = self.start_end(queue)[0]
                    qscpn = dict(zip(
                        qscpn.keys(), (
                            queue,
                            stack, 
                            children, 
                            pos, 
                            node,
                        )
                    ))
                    self.list_str_tree(s, **qscpn)
                else:
                    node.pop_child(children, qs)
                    child = node.add_child(qs, qe)
                    stack.append(queue.pop())
                    pos = self.start_end(queue)[0]
                    qscpn = dict(zip(
                        qscpn.keys(), (
                            queue, 
                            stack, 
                            children,
                            pos,
                            child,
                        )
                    ))
                    self.list_str_tree(s, **qscpn)
            elif qe > se:
                stack.pop()
                qscpn = dict(zip(
                    qscpn.keys(), (
                        queue, 
                        stack, 
                        children, 
                        pos, 
                        node.parent,
                    )
                ))
                self.list_str_tree(s, **qscpn)
        if children:
            if node.is_child(
                    qs, qe, ss, 
                    se
                ):
                if node.is_leaf(children, pos):
                    node.pop_child(children, qs)
                    node.add_leaf(qs, qe, s)
                    queue.pop()
            else:
                stack.pop()
                qscpn = dict(zip(
                    qscpn.keys(), ( 
                        queue, 
                        stack,
                        children, 
                        pos, 
                        node.parent
                    )
                ))
                self.list_str_tree(s, **qscpn)
        return node

    def list_str_tree_to_list(
            self, 
            tree, 
            lst,
            ):
        self.__list_str_tree_to_list(tree, lst)
        return lst

    
    def __list_str_tree_to_list(
            self, 
            node, 
            lst,
            ):
        if node.get_children():
            for child in node.get_children():
                if isinstance(child, type(node)):
                    self.utils.nested_list_insert(lst, [], node.path)
                    self.__list_str_tree_to_list(child, lst)
                if isinstance(child, list):
                    self.utils.nested_list_insert(lst, child, node.path)
        
    def str_to_list(
            self, 
            s,
            ):
        tree = self.list_str_tree(s)
        lst = []
        if self.utils.num_lists(s) == 1:
            return self.list_str_tree_to_list(tree, lst)[0]
        return self.list_str_tree_to_list(tree, lst)
    
def strcast(
        s,
        ):
    str_cast = StrCast()
    return str_cast.utils.str_cast(s)

