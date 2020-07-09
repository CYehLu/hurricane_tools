import os


###
### define some exception files
###

# the files end with '.py' and not start with '_' would be parsed into .md files
# if any exception file, list here
ignore_files = []
include_files = []


###
### define some functions
### 

def parse_one_line_doc(i, content):
    tri_quote = '"""'
    line = content[i]
    doc = line.strip().strip(tri_quote)
    return i+1, doc


def parse_multi_line_doc(i, content):
    tri_quote = '"""'
    i += 1
    doc_list = []
    while True:
        line = content[i]
        
        if line.strip() != tri_quote:
            doc_list.append(line)
            i += 1
        else:
            full_doc = ''.join(doc_list)
            
            strip_doc_list = list(map(lambda s: s.strip(), doc_list))
            blank_idx = strip_doc_list.index('')            
            func_descrpt = ' '.join(strip_doc_list[:blank_idx])
            break
            
    return i+1, func_descrpt, full_doc

    
    
def parse_pyfile(content):
    tri_quote = '"""'
    func_dict = {}   # key: function name, value: [function args, function description]
    doc_dict = {}    # key: function name, value: function full docstring

    i = 0
    in_func_scope = False
    while i < len(content):
        line = content[i]

        if line.startswith('def'):
            func_name = line.split()[1].split('(')[0]    # only function name
            func_name_args = line.strip('def').strip().strip(':')    # with args

            # move to next line
            i += 1
            in_func_scope = True
            continue

        if in_func_scope and tri_quote in line:

            # determine the docstring is one line or multiple lines
            # e.g: 
            # def oneline():
            #     """this is one line docstring"""
            #     pass
            # def multi_line():
            #    """
            #    this is multiple line docstring
            #    line2
            #    line3
            #    """
            #    pass

            if line.strip() == tri_quote:
                # parse multiline docstring
                i, func_descrpt, full_doc = parse_multi_line_doc(i, content)
                func_dict[func_name] = func_descrpt
                doc_dict[func_name] = [func_name_args, full_doc]
            else:
                # parse one line docstring
                i, doc = parse_one_line_doc(i, content)
                func_dict[func_name] = doc
                doc_dict[func_name] = [func_name_args, doc]

            in_func_scope = False

        elif in_func_scope and tri_quote not in line:
            # no docstring
            func_dict[func_name] = None
            doc_dict[func_name] = [func_name_args, None]
            in_func_scope = False
            i += 1

        if not in_func_scope:
            i += 1
            
    return func_dict, doc_dict


###
### start to parse .py files
###

# key: .py filename, value: dictionary with key is function name
module_dict = {
    'description': {},
    'fulldoc': {}
}
# e.g
# module_dict = {
#     'description': {
#         'file_a.py': {
#             'func_a1': 'The description of func_a1',
#             'func_a2': 'The description of func_a2'
#         },
#         'file_b.py': {
#             'func_b1': 'The description of func_b1'
#         }
#     },
#     'fulldoc': {
#         'file_a.py': {
#             'func_a1': ['func_a1(*args, **kwargs)', Full docstring of func_a1'],
#             'func_a2': ['func_a2(*args, **kwargs)', Full docstring of func_a2']
#         },
#         'file_b.py': {
#             'func_b1': ['func_b1(*args, **kwargs)', Full docstring of func_b1']
#         }
#     }
# }

pyfilenames = [f for f in os.listdir('./') if f.endswith('py') and not f.startswith('_')]
pyfilenames = pyfilenames + include_files
pyfilenames = [pyf for pyf in pyfilenames if pyf not in ignore_files]
pyfilenames = sorted(pyfilenames)

for pyfname in pyfilenames:
    with open(pyfname, 'r') as pyf_content:
        content = pyf_content.readlines()
        func_dict, doc_dict = parse_pyfile(content)
        module_dict['description'][pyfname] = func_dict
        module_dict['fulldoc'][pyfname] = doc_dict
        

###
### write parsed result into README.md
###  
    
with open('README.md', 'w') as mdfile:
    write_content = [
        '# hurricane_tools',
        'documents of all modules -> `./documents`',
        '\n'
    ]
    
    for module_name, sub_dict in module_dict['description'].items():
        write_content.append(module_name)
        write_content.append('------')
        write_content.append(f'[document](./documents/{module_name.rstrip(".py")}.md) \n')
        write_content.append('| Function | Description |')
        write_content.append('| :------- | :---------- |')
        
        for func_name, func_doc in sub_dict.items():
            write_content.append(f'| <font color="#a77864"> **{func_name}** </font> | {func_doc} |')
        
        write_content.append('\n')
        write_content.append('******')
        
    write_string = '\n'.join(write_content)
    mdfile.write(write_string)
    
    
###
### write parsed result into documents
### 

module_dict['description']
for pyfname in pyfilenames:
    with open(f'./documents/{pyfname.rstrip(".py")}.md', 'w') as mdfile:
        write_content = [f'# {pyfname}']
        
        for func_name_args, doc in module_dict['fulldoc'][pyfname].values():
            # e.g func_name_args = 'func1(arg1, arg2, arg3)'
            name = func_name_args.split('(')[0]
            args = '(' + func_name_args.split('(')[1]
            write_content.append(f'<span style="color:#a77864">**{name}**</span>**{args}**\n')
            write_content.append(doc)
            write_content.append('\n')
            write_content.append('******')
            
        write_string = '\n'.join(write_content)
        mdfile.write(write_string)

