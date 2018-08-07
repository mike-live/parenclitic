import timeit
import numpy as np
import os.path
import sys
import os
from pathlib2 import Path
from configparser import ConfigParser

class param:
    def __init__(self, value_be = None, value_en = None, num_ticks = 1, id_tick = 0, whole_values = False, manual_ticks = False, value_list = None, dtype = None, name = '', is_path = True):
        self.value_be = value_be
        if value_en is None:
            self.value_en = value_be
        else:    
            self.value_en = value_en
        self.num_ticks = num_ticks
        if dtype is None:
            self.dtype = type(value_be)
        else:
            self.dtype = dtype
        self.is_const = num_ticks == 1
        self.id_tick = id_tick
        self.manual_ticks = manual_ticks
        self.whole_values = whole_values
        if value_be is None and value_list is not None:
            self.value = value_list[0];
        else:
            self.value = value_be
        self.name = name
        self.is_path = is_path        
        self.value_list = value_list

    def set_tick(self, id_tick):
        self.id_tick = id_tick
        if id_tick >= self.num_ticks or id_tick < 0:
            raise IndexError
        if isinstance(self.value_be, (np.floating, float, int, long)) and isinstance(self.value_en, (np.floating, float, int, long)):
            if self.num_ticks > 1: self.value = (self.value_en - self.value_be) * id_tick / (self.num_ticks - 1) + self.value_be
            else: self.value = value_be
        elif self.value_list and hasattr(self.value_list, '__getitem__'):
            self.value = self.value_list[id_tick]
        return self.value
            
    def set_value(self, value):
        self.value = value

    def get_values(self):
        values = map (self.set_tick, range(self.num_ticks))
        return values

    def __iter__(self):
        self.id_tick = -1
        return self

    def next(self):
        try:
            value = self.set_tick(self.id_tick + 1)
        except IndexError:
            raise StopIteration
        return value

    def __len__(self):
        return self.num_ticks

    def __getitem__(self, id_tick = 0):
        return self.set_tick(id_tick)

def set_params(params, id_tick):
    for name, param in params.iteritems():
        if not param.is_const and not param.manual_ticks:
            param.set_tick(id_tick % param.num_ticks)
            id_tick /= param.num_ticks
    

class configuration:
    def __init__(self, params, info, files, data_path = 'data', data_name = '', project_name = '', config_name = '', params_sets = {}):
        self.params = params
        #self.init_type = init_type
        self.info = info
        self.files = files
        self.run_id = info['run_id']
        self.run_num = info['run_num']
        self.project_path = info['work_dir']
        for parent in Path(info['work_dir']).parents:
            if parent.stem == project_name:
                self.project_path = parent
                break
        self.upd_ticks()
        self.data_path = self.project_path / data_path / data_name
        self.name = config_name
        self.params_sets = params_sets

    def upd_ticks(self):
        set_params(self.params, self.run_id)

    def get_joined_params(self, params_str = '', include_set = None, exclude_set = None, need_const = False, need_nonconst = False, delimiter = '_'):
        for param_name, param in self.params.iteritems():
            need_param = ((param.is_const and need_const) or (not param.is_const and need_nonconst)) and param.is_path
            if include_set:
                need_param &= param_name in include_set
            if exclude_set:
                need_param &= param_name not in exclude_set
            #print param_name, need_param, param.is_const, param.is_path, need_const, need_nonconst
            if need_param:
                if isinstance(params_str, Path):
                    params_str /= param_name + delimiter + str(param.value)
                else:
                    if params_str:
                        params_str += delimiter
                    params_str += param_name + delimiter + str(param.value)
        #print need_nonconst, params_str
        return params_str

    def get_path_params(self, name, include_set = None, exclude_set = None, need_const = True, need_nonconst = True, delimiter = '_'):
        dir_name = self.get_joined_params(Path(), include_set, exclude_set, need_const = need_const, delimiter = delimiter)
        name = Path(name)
        name = name.with_name(self.get_joined_params(name.stem, include_set, exclude_set, need_nonconst = need_nonconst, delimiter = delimiter) + name.suffix)
        path_const_params = self.data_path / 'params' / dir_name / name
        if not path_const_params.parent.exists():
            path_const_params.parent.mkdir(parents = True, exist_ok = True)
        return path_const_params

    def ofname(self, names, ext = '', include_set = None, exclude_set = None, delimiter = '_'):
        path = Path('')
        if type(names) is list:
            for name in names:
                if type(name) is list:
                    path /= delimiter.join(map(lambda xname: str(self.files.get(xname, xname)), name))
                else:
                    path /= self.files.get(name, name)
        else:
            path /= self.files.get(names, names)
        return str(self.get_path_params(path.with_name(path.name + ext), include_set, exclude_set, delimiter = delimiter))

    def ifname(self, name):
        cur_path = self.data_path / self.files[name]
        return str(cur_path.resolve())

    def save_params(self, path = None, include_set = None, exclude_set = None):
        if path is None:
            path = self.get_path_params('config.ini', include_set, exclude_set, need_nonconst = False)

        config_struct = ConfigParser()
        
        section_name = 'Run information'
        config_struct.add_section(section_name)
        for name, value in self.info.iteritems():
            config_struct.set(section_name, name, str(value))

        section_name = 'Parameters'
        config_struct.add_section(section_name)
        for name, param in self.params.iteritems():
            if param.num_ticks > 1:
                if param.value_list is None:
                    param_values = "linspace(" + str(param.value_be) + ", " + str(param.value_en) + ", " + str(param.num_ticks) + ")"
                else:
                    param_values = str(param.value_list)

            if param.whole_values:
                param_str = param_values
            else:
                param_str = str(param.value)
                if not param.is_const:
                    param_str += " in " + param_values

            config_struct.set(section_name, param.name, param_str)

        section_name = 'Files'
        config_struct.add_section(section_name)
        for name, fname in self.files.iteritems():
            config_struct.set(section_name, name, str(fname))
        
        config_struct.write(Path(path).open("w"))
