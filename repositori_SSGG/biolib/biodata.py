#!/usr/bin/env python
"""
SYNOPSIS
    biodata.py [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    biodata.py: Module including abstract class definitions for any biological data, sets of biological data, rendering classes, etc.
    All future biological data classes should be derived from this one.
EXAMPLES
    biodata.py
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    biodata.py v0.1 (from September 2013)
"""

from abc import ABCMeta, abstractmethod, abstractproperty
import collections
import sys
from overloading import Overloaded


class BioContainer(collections.MutableSet):
    def __init__(self, iterable=None):
        self.end = end = []
        end += [None, end, end]         # sentinel node for doubly linked list
        self.map = {}                   # key --> [key, prev, next]
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def add(self, key):
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]

    def discard(self, key):
        if key in self.map:
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = self.end[1][0] if last else self.end[2][0]
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, BioContainer):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)


class BioData:
    __metaclass__ = ABCMeta
    pass


class Variant(BioData):
    """Base class for storing variant information"""
    def __init__(self):
        self._position_info = {}
        self._samples_info = {}
        self._annotations_info = {}

    @Overloaded
    def __repr__(self, render_obj):
        """TODO: find a way to define the repr method, so that it gets an after renderer sorted list of parameters"""
        render_obj.render(self)

    def __iter__(self):
        return self

    @Overloaded
    def get_data(self, render_obj, input_obj):
        render_obj.get_data(self, input_obj)

    @property
    def position_info(self):
        """Position related information for variant: chromosome, position, reference, variant"""
        return self._position_info

    @position_info.setter
    def position_info(self, position_info):
        self._position_info = position_info

    @position_info.deleter
    def position_info(self):
        del self._position_info

    @property
    def samples_info(self):
        """Samples related information for variant: depth, allele frequency, zygosity or genotype"""
        return self._samples_info

    @samples_info.setter
    def samples_info(self, samples_info):
        self._samples_info = samples_info

    @samples_info.deleter
    def samples_info(self):
        del self._samples_info

    @property
    def annotations_info(self):
        """Annotations related information for variant: dependent on deployed databases and versions"""
        return self._annotations_info

    @annotations_info.setter
    def annotations_info(self, annotations_info):
        self._annotations_info = annotations_info

    @annotations_info.deleter
    def annotations_info(self):
        del self._annotations_info

BioData.register(Variant)


class Renderer:
    """esta clase deberia ser generica, incluyendo la clase a representar y otros atributos"""
    __metaclass__ = ABCMeta
    pass


class SeparatedValuesRenderer(Renderer):
    def __init__(self, separator='|', stream=None, indent=0, **kwargs):
        self.separator = kwargs.get('header') if kwargs.get('header') is not None else separator
        self.stream = stream if stream is not None else sys.stdout
        self.indent = indent
        self.config = []

    def write(self, text):
        for line in text.splitlines(True):
            self.stream.write(line)

    def get_template(self, file_path='', **kwargs):
        """Two types of construction: either gets configuration of a results psv header
        or values are directly assigned"""

        # Check tuples and keyed-arguments. If "header" argument is found, create from line:
        file_path = kwargs.get('file_path') if kwargs.get('file_path') is not None else file_path
        file_handle = open(file_path, 'r')

        for line in file_handle.readlines():
            if line.startswith('#'):
                header = line[1:]
                tmp = header.split(self.separator)
                current_sample = ''
                for field in tmp:
                    if field.endswith('_Genotype') or field.endswith('_Depth') or field.endswith('_Var/Depth'):
                        sample_field = ''
                        if field.endswith('_Genotype'):
                            sample_name = field[0:-9]
                            sample_field = 'Genotype'
                        elif field.endswith('_Depth'):
                            sample_name = field[0:-6]
                            sample_field = 'Depth'
                        elif field.endswith('_Var/Depth'):
                            sample_name = field[0:-10]
                            sample_field = 'Var/Depth'
                        if current_sample == '' or current_sample != sample_name:
                            current_sample = sample_name
                            if sample_name == 'position_info' or sample_name == 'annotations_info':
                                print 'ERROR: Invalid sample name',sample_name
                                raise
                        self.config.append((sample_name, sample_field))
                    else:
                        if current_sample == '':
                            self.config.append(('position_info', field))
                        else:
                            self.config.append(('annotations_info', field))
            else:
                file_handle.close()
                return

    @Overloaded
    def render(self, obj):
        text = ''
        for parameter in self.config:
            if hasattr(obj, parameter[0]):
                attribute = getattr(obj, parameter[0])
                text += attribute[parameter[1]]
            else:
                attribute = getattr(obj, 'samples_info')
                text += attribute[parameter[0]][parameter[1]]
            text += self.separator
        text = text[:-1]
        self.write(text)
        return text

    @render.register(object, list)
    def render_list(self, obj):
        if not obj:
            self.write("[]")
            return
        sep = "["
        self.indent += 1
        for item in obj:
            self.write(sep)
            self.render(item)
            sep = ",\n"
        self.indent -= 1
        self.write("]")

    @Overloaded
    def get_data(self, output_obj, input_obj):
        tmp = input_obj.split(self.separator)
        column = 0
        for field in self.config:
            if field[0] == 'position_info':
                output_obj.position_info[field[1]] = tmp[column]
            elif field[0] == 'annotations_info':
                output_obj.annotations_info[field[1]] = tmp[column]
            else:
                if field[0] in output_obj.samples_info:
                    output_obj.samples_info[field[0]][field[1]] = tmp[column]
                else:
                    output_obj.samples_info[field[0]] = {field[1]: tmp[column]}
            column += 1


class HexPrettyPrinter(SeparatedValuesRenderer):
    @Overloaded
    def render(self, obj):
        SeparatedValuesRenderer.render(self, obj)

    @render.register(object, int)
    @render.register(object, long)
    def pprint_hex(self, obj):
        self.write(hex(obj))


class Rendering:
    """esta clase deberia ser generica, incluyendo la clase a representar y otros atributos"""
    __metaclass__ = ABCMeta
    pass


class SeparatedValuesRendering(Rendering):
    """Base class for rendering BioData and BioContainers into/from lists of separated values"""
    def __init__(self, *args, **kwargs):
        # Check tuples and keyed-arguments:
        self._separator = args[0] if args else None
        if self._separator is None:
            self._separator = kwargs.get('header')

    @property
    def separator(self):
        return self._separator

    @separator.setter
    def separator(self, separator):
        self._separator = separator

    def get_header(self, variant, *args, **kwargs):
        """Two types of construction: either gets configuration of a results psv header
        or values are directly assigned"""

        # Check tuples and keyed-arguments. If "header" argument is found, create from line:
        header = args[0] if args else None
        if header is None:
            header = kwargs.get('header')

        tmp = header.split(self.separator)
        current_sample = ""
        for field in tmp:
            if field.endswith('_Genotype') or field.endswith('_Depth') or field.endswith('_Var/Depth'):
                sample_field = ""
                if field.endswith('_Genotype'):
                    sample_name = field[0:-9]
                    sample_field = 'Genotype'
                elif field.endswith('_Depth'):
                    sample_name = field[0:-6]
                    sample_field = 'Depth'
                elif field.endswith('_Var/Depth'):
                    sample_name = field[0:-10]
                    sample_field = 'Var/Depth'
                if current_sample == "" or current_sample != sample_name:
                    current_sample = sample_name
                    variant.samples_info.update({sample_name: {}})
                variant.samples_info[sample_name].update({sample_field: ''})
            else:
                if current_sample == '':
                    variant.position_info.update({field: ''})
                else:
                    variant.annotations_info.update({field: ''})

    def get_variant_info(self, variant, *args, **kwargs):
        """Two types of construction: either gets configuration of a results psv header
        or values are directly assigned"""

        # Check tuples and keyed-arguments. If "header" argument is found, create from line:
        line = args[0] if args else None
        if line is None:
            line = kwargs.get('header')

        values = line.split(self.separator)
        mycol = 0
        for index, field in enumerate(self.variant_info):
            field_name = self.variant_info[index].keys()[0]
            self.variant_info[index][field_name] = values[mycol]
            mycol += 1
        for index, sample in enumerate(self.samples_info):
            sample_name = self.samples_info[index].keys()[0]
            for samp_index, field in enumerate(self.samples_info[index][sample_name]):
                field_name = self.samples_info[index][sample_name][samp_index].keys()[0]
                self.samples_info[index][sample_name][samp_index][field_name] = values[mycol]
                mycol += 1
        for index, field in enumerate(self.annotations_info):
            field_name = self.annotations_info[index].keys()[0]
            self.annotations_info[index][field_name] = values[mycol]
            mycol += 1

SeparatedValuesRendering.register(BioData)
SeparatedValuesRendering.register(BioContainer)

class FutureIterator:
    def next(self):
        self.__state__ += 1
        if self.__state__ > len(self.position_info) + len(self.samples_info) + len(self.annotations_info):
            raise StopIteration
        elif self.__state__ > len(self.position_info) + len(self.samples_info):
            return
        if self.current > self.high:
            raise StopIteration
        else:
            self.current += 1
            return self.current - 1