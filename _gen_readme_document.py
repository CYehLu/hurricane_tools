import os


def parse_one_line_doc(i, content):
    """
    parse one line docstring of a function
    
    input: i and content
        content[i] = 'def func(args):\n'
        content[i+1] = '\"""this line is one-line docstring\n\"""'
        content[i+2] = '# this line and below is function body'
    output: i+2 and docstring
    """
    tri_quote = '"""'
    
    # adding 4 blank in order to be consistent with multi-line situation
    doc = '    ' + content[i+1].strip().strip(tri_quote)
    return i+2, doc 


def parse_multi_line_doc(i, content):
    """
    parse multi line docstring of a function
    
    input: i and content
        content[i] = 'def func(args):\n'
        content[i+1] = '\"""\n'   (triple quote)
        content[i+2] = 'function description line 1\n'  -| (function description)  -|
        content[i+3] = 'function description line 2\n'  -|                          |
        content[i+4] = '\n'                                                         | (full doc)
        contnet[i+5] = 'still docstring\n'              -|                          |
        ...                                              | (addition part)          |
        content[i+n] = 'last line of docstring\n'       -|                         -|
        content[i+n+1] = '\"""\n'
        content[i+n+2] = '# function body start\n'
    output: i+n+2, function description, and full docstring
    """
    tri_quote = '"""'
    i = i + 2
    
    func_descrpt = []
    full_doc = []
    
    # here is function description part   
    while content[i].strip() and content[i].strip() != tri_quote:
        func_descrpt.append(content[i].strip())
        full_doc.append(content[i])
        i += 1
        
    # check if has blank line, or has reached the triple quote
    if not content[i].strip():
        full_doc.append(content[i])
        i += 1
        
        # addition part (parameters, retruns, note, reference ... etc)
        while content[i].strip() != tri_quote:
            full_doc.append(content[i])
            i += 1

        # now `i` is at the end triple quote of docstring
        # convert `func_descrpt` and `full_doc` from list to string
        func_descrpt = ' '.join(func_descrpt)
        full_doc = ''.join(full_doc)

        return i+1, func_descrpt, full_doc
    
    else:
        # reach triple quote -> docstring only has function description part
        func_descrpt = ' '.join(func_descrpt)
        full_doc = ''.join(full_doc)
        return i+1, func_descrpt, full_doc
    

def parse_func(i, content):
    """
    content[i] = 'def func_name(args):\n'
    return 1. function body index
           2. function name
           3. name+arguments
           4. function description
           5. full doc
    example: (
        10,
        'func_name', 
        'func_name(arg1, arg2)', 
        'this is func description line1\nline 2',
        'this is func description line1\nline 2\n\nparameters:...\n'
    )
    """
    tri_quote = '"""'
    
    # parse name and arg
    func_name = content[i].strip().split()[1].split('(')[0]
    func_name_args = content[i].strip().split('def')[1].split(':')[0].strip()
    
    # determine one-line or multi-line docstring
    if tri_quote not in content[i+1]:
        # there is no docstring
        func_descrpt = ''
        full_doc = ''
        i += 1
        
    elif content[i+1].strip() == tri_quote:
        i, func_descrpt, full_doc = parse_multi_line_doc(i, content)
        
    else:
        i, doc = parse_one_line_doc(i, content)
        func_descrpt = doc
        full_doc = doc
    
    return i, func_name, func_name_args, func_descrpt, full_doc
    

def parse_cls(i, content):
    """
    return: tuple
        (1) index of the first line of non-class part
        (2) class name, 
        (3) class description, 
        (4) detail class description (cls_full_doc), 
        (5) method_dict
            and method_dict = {
                method_name: [name_and_arg, method_description, method_full_doc]
            }
    example
        [12] class Car:
        [13]     \"""
        [14]     Car description
        [15]
        [16]     Some detail information of this class
        [17]     ...still the detail information
        [18]     \"""
        [19]     def __init__(self, arg1, arg2):
        [20]         \"""one-line docstring for __init__\"""
        [21]         pass
        [22]
        [23]
        [24] def a_function(args):
        [25]     pass
        
        then parse_cls would return
        (24, 
         'Car', 
         'Car description', 
         'Car description\n\nSome detail information...',
         {'__init__': ('__init__(self, arg1, arg2)',
                       'one-line docstring for __init__',
                       'one-line docstring for __init__')
         }
        )
    """
    tri_quote = '"""'
    
    # parse class name and description
    cls_name = content[i].split()[1].split('(')[0].split(':')[0]
    
    if content[i+1].strip().startswith(tri_quote):
        # there is a docstring for class
        # use the same way of parse function docstring to get the class description
        if content[i+1].strip() == tri_quote:
            i, cls_descrpt, cls_full_doc = parse_multi_line_doc(i, content)
        else:
            i, cls_descrpt = parse_one_line_doc(i, content)
            cls_full_doc = cls_descrpt
            
    else:
        # there is no docstring for this class
        cls_descrpt = ''
        cls_full_doc = ''
        i += 1
        
    # parse methods of this class
    # do not parse private methods start with '_' or '__', apart from '__init__' and '__call__'
    allowed_underscore = ['__init__', '__call__']
    method_dict = {}
    in_cls_scope = True
    
    while i < len(content):           
        in_cls_scope = content[i] and (not content[i][0].isalpha()) and (not content[i][0].isdigit())
        
        if in_cls_scope:
            # still in the class scope
            
            if content[i].strip().startswith('def'):
                res = parse_func(i, content)
                i, method_name, method_name_args, method_descrpt, method_full_doc = res
                
                if method_name.startswith('_') and method_name not in allowed_underscore:
                    # not to parse private methods
                    continue
                else:
                    method_dict[method_name] = (
                        method_name_args,
                        method_descrpt,
                        method_full_doc
                    )
                    
            else:
                i += 1
            
        else:
            # have leaved class scope
            break
            
    return i, cls_name, cls_descrpt, cls_full_doc, method_dict


def parse_module(file_path):
    module_name = file_path.split('/')[-1].rstrip('.py')
    module_dict = {}   # key is function/class name, value is some informations
    
    ### e.g
    ### module_dict = {
    ###     'func1': [
    ###         '<function>'
    ###         'func1(arg1, arg2)', 
    ###         'func1_description', 
    ###         'func1_description\nParameters:...(full_doc)'
    ###      ]
    ###     'func2': [
    ###         '<function>'
    ###         'func2(arg1, arg2)', 
    ###         'func2_description', 
    ###         'func2_description\nParameters:...(full_doc)'
    ###      ]
    ###     'cls1': (
    ###         '<class>'
    ###         'class description',
    ###         'more detail class description (similar to `full_doc` of function)'
    ###         {'__init__': [
    ###              '__init__(self, arg1, arg2)',
    ###              '__init__description',
    ###              '__init__description\nParameters:...(full_doc)'
    ###          ],
    ###          'method1': [
    ###              'method1(self, arg1, arg2)',
    ###              'method1_description',
    ###              'method1_description\nParamters:...(full_doc)'
    ###          ],
    ###          ...
    ###         }
    ###     )
    ### }
    
    with open(file_path) as pyf:
        content = pyf.readlines()
        
    i = 0
    while i < len(content):
        line = content[i]
        
        if line.startswith('def'):
            if line.split()[1].startswith('_'):
                # this is private function, e.g: 'def _private_func():\n'
                # not to parse private function
                i += 1
                continue
            else:
                res = parse_func(i, content)
                i, func_name, func_name_args, func_descrpt, full_doc = res
                module_dict[func_name] = [
                    '<function>', 
                    func_name_args, 
                    func_descrpt, 
                    full_doc
                ]
                
        elif line.startswith('class'):
            if line.split()[1].startswith('_'):
                # private class
                i += 1
                continue
            else:
                res = parse_cls(i, content)
                i, cls_name, cls_descrpt, cls_full_doc, method_dict = res
                module_dict[cls_name] = (
                    '<class>',
                    cls_descrpt, 
                    cls_full_doc, 
                    method_dict
                )
        
        else:
            i += 1
            
    return module_name, module_dict


def parse_pkgs(pkgs_path):
    """return dict"""
    pass


class ReadmeContent:
    def __init__(self):
        self.content = [
            '# hurricane_tools',
            'documents of all modules -> `./documents`',
            '\n'
        ]
        
    def update_content(self, module_name, module_dict):
        self.content.append(module_name)
        self.content.append('------')
        self.content.append(f'[document](./documents/{module_name}.md) \n')
        self.content.append('| Function / Class | Description |')
        self.content.append('| :--------------- | :---------- |')
        
        for func_cls_name, container in module_dict.items():
            if container[0] == '<function>':
                description = container[2]
            elif container[0] == '<class>':
                description = container[1]
                
            self.content.append(f'| <font color="#a77864"> **{func_cls_name}** </font> | {description.strip()} |')
        
        self.content.append('\n')    
        self.content.append('******')
            
    def write(self):
        with open('README.md', 'w') as file:
            write_content = '\n'.join(self.content)
            file.write(write_content)


class DocumentContent:
    @staticmethod
    def _append_function(name, container, write_content):
        name_args = container[1]
        args = '(' + name_args.split('(')[1]
        full_doc = container[3]
        
        if '**' in args:
            args = args.replace('*', '\*')
        
        write_content.append(f'<span style="color:#a77864">**{name}**</span>**{args}**\n')
        write_content.append(full_doc)
        write_content.append('\n')
        write_content.append('******')
    
    @staticmethod
    def _append_class(name, container, write_content):
        # class description
        full_doc = container[2]
        write_content.append(f'class <span style="color:#a77864">**{name}**</span>\n')
        write_content.append(full_doc)
        write_content.append('\n')
                
        # create a tabel to list all class methods and their description
        method_dict = container[3]
        
        write_content.append('| Methods | Description |')
        write_content.append('| :------ | :---------- |')
        
        for method_name, method_info in method_dict.items():
            method_name = method_name.replace('_', '\_')
            write_content.append(f'| <font color="#a77864"> **{method_name}** </font> | {method_info[1].strip()} |')
        write_content.append('\n')
        
        # detail of methods
        method_dict = container[3]
        for method_name, method_info in method_dict.items():
            name_args = method_info[0]
            args = '(' + name_args.split('(')[1]
            full_doc = method_info[2]
            
            # insert '\' before every underscore '_'
            #method_name = ''.join(['\{}'.format(c) if c == '_' else c for c in method_name])
            method_name = method_name.replace('_', '\_')
            
            write_content.append(f'<span style="color:#cca99b">{name}</span>.<span style="color:#a77864">**{method_name}**</span>**{args}**\n')
            write_content.append(f'{full_doc}')
            write_content.append('  ')
        
        write_content.append('******')
    
    @classmethod
    def write_document(cls, module_name, module_dict):
        write_content = [
            f'# {module_name}  \n',
            f'[[source](../{module_name}.py)]  \n'
        ]

        for func_cls_name, container in module_dict.items():
            if container[0] == '<function>':
                cls._append_function(func_cls_name, container, write_content)
            elif container[0] == '<class>':
                cls._append_class(func_cls_name, container, write_content)
                
        write_content = '\n'.join(write_content)

        if not os.path.isdir('./documents'):
            os.mkdir('./documents')
            
        with open(f'./documents/{module_name}.md', 'w') as file:
            file.write(write_content)
            

def gen_readme_document(include=None, ignore=None):
    if include is None:
        include = []
    if ignore is None:
        ignore = []
        
    ls = os.listdir()
    
    is_folder = [True if os.path.isdir(l) else False for l in ls]
    is_py = [True if l.endswith('.py') else False for l in ls]
    is_private = [True if l.startswith('_') else False for l in ls]
    is_ignore = [True if l in ignore else False for l in ls]
    
    cond = [is_folder[i] or is_py[i] and not is_private[i] and not is_ignore[i] for i in range(len(ls))]
    ls = include + [ls[i] for i in range(len(ls)) if cond[i]]
    ls = sorted(ls)
    
    readme = ReadmeContent()
    
    for ff in ls:
        if os.path.isdir(ff):
            parse_pkgs(None)
        else:
            module_name, module_dict = parse_module('./' + ff)
            readme.update_content(module_name, module_dict)
            DocumentContent.write_document(module_name, module_dict)
            
    readme.write()


if __name__ == '__main__':
    include = []
    ignore = []
    gen_readme_document(include, ignore)