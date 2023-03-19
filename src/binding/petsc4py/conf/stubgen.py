import inspect
from textwrap import dedent

def is_cyfunction(obj):
    return type(obj).__name__ == 'cython_function_or_method'


def is_function(obj):
    return (
        inspect.isbuiltin(obj)
        or is_cyfunction(obj)
        or type(obj) is type(ord)
    )


def is_method(obj):
    return (
        inspect.ismethoddescriptor(obj)
        or inspect.ismethod(obj)
        or is_cyfunction(obj)
        or type(obj) in (
            type(str.index),
            type(str.__add__),
            type(str.__new__),
        )
    )


def is_classmethod(obj):
    return (
        inspect.isbuiltin(obj)
        or type(obj).__name__ in (
            'classmethod',
            'classmethod_descriptor',
        )
    )


def is_staticmethod(obj):
    return (
        type(obj).__name__ in (
            'staticmethod',
        )
    )

def is_constant(obj):
    return isinstance(obj, (int, float, str))

def is_datadescr(obj):
    return inspect.isdatadescriptor(obj) and not hasattr(obj, 'fget')


def is_property(obj):
    return inspect.isdatadescriptor(obj) and hasattr(obj, 'fget')


def is_class(obj):
    return inspect.isclass(obj) or type(obj) is type(int)


class Lines(list):

    INDENT = " " * 4
    level = 0

    @property
    def add(self):
        return self

    @add.setter
    def add(self, lines):
        if lines is None:
            return
        if isinstance(lines, str):
            lines = dedent(lines).strip().split("\n")
        indent = self.INDENT * self.level
        for line in lines:
            self.append(indent + line)


def signature(obj):
    doc = obj.__doc__
    doc = doc or f"{obj.__name__}: Any"  # FIXME remove line
    sig = doc.split('\n', 1)[0].split('.', 1)[-1]
    return sig

def stubgen_constant(constant):
    name, value = constant
    return f"{name}: Final[{type(value).__name__}] = ..."

def stubgen_function(function):
    sig = signature(function)
    return f"def {sig}: ..."

def stubgen_method(method):
    sig = signature(method)
    return f"def {sig}: ..."

def stubgen_datadescr(datadescr):
    sig = signature(datadescr)
    return f"{sig}"

def stubgen_property(prop, name=None):
    sig = signature(prop.fget)
    pname = name or prop.fget.__name__
    ptype = sig.rsplit('->', 1)[-1].strip()
    return f"{pname}: {ptype}"

def stubgen_constructor(cls, name='__init__'):
    init = (name == '__init__')
    argname = cls.__name__.lower()
    argtype = cls.__name__
    initarg = f"{argname}: Optional[{argtype}] = None"
    selfarg = 'self' if init else 'cls'
    rettype = 'None' if init else argtype
    arglist = f"{selfarg}, {initarg}"
    sig = f"{name}({arglist}) -> {rettype}"
    return f"def {sig}: ..."

def stubgen_class(cls, outer=None, done=None):
    skip = {
        '__doc__',
        '__dict__',
        '__module__',
        '__weakref__',
        '__pyx_vtable__',
        '__enum2str',  # FIXME refactor implemetation
        '_traceback_', # FIXME maybe refactor?
        '__lt__',
        '__le__',
        '__ge__',
        '__gt__',
    }
    special = {
        '__len__':   "__len__(self) -> int",
        '__bool__':  "__bool__(self) -> bool",
        '__hash__':  "__hash__(self) -> int",
        '__int__':   "__int__(self) -> int",
        '__index__': "__int__(self) -> int",
        '__str__':   "__str__(self) -> str",
        '__repr__':  "__repr__(self) -> str",
        '__eq__':    "__eq__(self, other: object) -> bool",
        '__ne__':    "__ne__(self, other: object) -> bool",
    }

    qualname = cls.__name__
    cls_name = cls.__name__
    if outer is not None and cls_name.startswith(outer):
        cls_name = cls_name[len(outer):]
        qualname = f"{outer}.{cls_name}"

    override = OVERRIDE.get(qualname, {})
    done = set() if done is None else done
    lines = Lines()

    try:
        class sub(cls):
            pass
        final = False
    except TypeError:
        final = True
    if final:
        lines.add = "@final"
    base = cls.__base__
    if base is object:
        lines.add = f"class {cls_name}:"
    else:
        lines.add = f"class {cls_name}({base.__name__}):"
    lines.level += 1
    start = len(lines)

    for name in ('__new__', '__init__', '__hash__'):
        if name in cls.__dict__:
            done.add(name)

    dct = cls.__dict__
    keys = list(dct.keys())

    def dunder(name):
        return name.startswith('__') and name.endswith('__')

    def members(seq):
        for name in seq:
            if name in skip:
                continue
            if name in done:
                continue
            if dunder(name):
                if name not in special and name not in override:
                    done.add(name)
                    continue
            yield name

    for name in members(keys):
        attr = getattr(cls, name)
        if is_class(attr):
            done.add(name)
            lines.add = stubgen_class(attr, outer=cls_name)
            continue

    for name in members(keys):

        if name in override:
            done.add(name)
            lines.add = override[name]
            continue

        if name in special:
            done.add(name)
            sig = special[name]
            lines.add = f"def {sig}: ..."
            continue

        attr = getattr(cls, name)

        if is_method(attr):
            done.add(name)
            if name == attr.__name__:
                obj = dct[name]
                if is_classmethod(obj):
                    lines.add = "@classmethod"
                elif is_staticmethod(obj):
                    lines.add = "@staticmethod"
                lines.add = stubgen_method(attr)
            else:
                lines.add = f"{name} = {attr.__name__}"
            continue

        if is_datadescr(attr):
            done.add(name)
            lines.add = stubgen_datadescr(attr)
            continue

        if is_property(attr):
            done.add(name)
            lines.add = stubgen_property(attr, name)
            continue

        if is_constant(attr):
            done.add(name)
            lines.add = stubgen_constant((name, attr))
            continue

    leftovers = [name for name in keys if
                 name not in done and name not in skip]
    if leftovers:
        raise RuntimeError(f"leftovers: {leftovers}")

    if len(lines) == start:
        lines.add = "pass"
    lines.level -= 1
    return lines


def stubgen_module(module, done=None):
    skip = {
        '__doc__',
        '__name__',
        '__loader__',
        '__spec__',
        '__file__',
        '__package__',
        '__builtins__',
        '__pyx_unpickle_Enum',  # FIXME review
    }

    done = set() if done is None else done
    lines = Lines()

    keys = list(module.__dict__.keys())
    keys.sort(key=lambda name: name.startswith("_"))

    constants = [
        (name, getattr(module, name)) for name in keys
        if (
            name not in done and name not in skip and
            isinstance(getattr(module, name), int)
        )
    ]
    for attr in constants:
        name, value = attr
        done.add(name)
        if name in OVERRIDE:
            lines.add = OVERRIDE[name]
        else:
            lines.add = stubgen_constant((name, value))
    if constants:
        lines.add = ""

    for name in keys:
        if name in done or name in skip:
            continue
        value = getattr(module, name)

        if is_class(value):
            done.add(name)
            if value.__module__ != module.__name__:
                continue
            lines.add = stubgen_class(value)
            lines.add = ""
            instances = [
                (k, getattr(module, k)) for k in keys
                if (
                    k not in done and k not in skip and
                    type(getattr(module, k)) is value
                )
            ]
            for attrname, attrvalue in instances:
                done.add(attrname)
                lines.add = stubgen_constant((attrname, attrvalue))
            if instances:
                lines.add = ""
            continue

        if is_function(value):
            done.add(name)
            if name == value.__name__:
                lines.add = stubgen_function(value)
            else:
                lines.add = f"{name} = {value.__name__}"
            continue

    lines.add = ""
    for name in keys:
        if name in done or name in skip:
            continue
        value = getattr(module, name)
        done.add(name)
        if name in OVERRIDE:
            lines.add = OVERRIDE[name]
        else:
            lines.add = stubgen_constant((name, value))

    leftovers = [name for name in keys if
                 name not in done and name not in skip]
    if leftovers:
        raise RuntimeError(f"leftovers: {leftovers}")
    return lines


IMPORTS = """
from __future__ import annotations
import sys
from threading import Lock
from typing import (
    Any,
    Union,
    Optional,
    NoReturn,
    overload,
)
if sys.version_info >= (3, 8):
    from typing import (
        final,
        Final,
        Literal,
    )
else:
    from typing_extensions import (
        final,
        Final,
        Literal,
    )
if sys.version_info >= (3, 9):
    from collections.abc import (
        Callable,
        Hashable,
        Iterable,
        Iterator,
        Sequence,
        Mapping,
    )
else:
    from typing import (
        Callable,
        Hashable,
        Iterable,
        Iterator,
        Sequence,
        Mapping,
    )
if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self
from os import PathLike

import numpy

IntType: numpy.dtype = ...
RealType: numpy.dtype = ...
ComplexType: numpy.dtype = ...
ScalarType: numpy.dtype = ...
"""

OVERRIDE = {
    'Error': {
    },
    '__pyx_capi__': "__pyx_capi__: Final[Dict[str, Any]] = ...",
    '__type_registry__': "__pyx_capi__: Final[Dict[int, type]] = ...",
}

TYPING = """
"""

def stubgen_petsc4py_PETSc(done=None):
    from petsc4py import PETSc

    lines = Lines()

    lines.add = IMPORTS

    lines.add = ""
    lines.add = TYPING

    lines.add = ""
    lines.add = stubgen_module(PETSc)

    return lines


def main():
    import os
    output = os.path.join('src', 'petsc4py', 'PETSc.pyi')
    with open(output, 'w') as f:
        for line in stubgen_petsc4py_PETSc():
            print(line, file=f)


if __name__ == '__main__':
    main()
