import os


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
    
    # dealing two types
    #     def first_type():
    #         """
    #         This is the first line of docstring.
    #         second line
    #         third line
    #         
    #         Parameter:
    #         ...
    #         """
    #         pass
    #
    #     def second_type():
    #         """
    #         This is the first line of docstring.
    #         second line
    #         """
    #         pass
    
    doc_list = []
    while True:
        line = content[i]
        if line.strip() and line.strip() != tri_quote:
            doc_list.append(line.strip())
            i += 1
        elif line.strip() and line.strip() == tri_quote:
            # second type, the end of docstring
            doc = ' '.join(doc_list)
            break
        elif not line.strip():
            # the blank line in the first type
            doc = ' '.join(doc_list)
            break
            
    return i+1, doc


def parse_pyfile(content):
    func_dict = {}    # key: function name, value: function doc-string

    tri_quote = '"""'
    in_func_scope = False
    i = 0

    while i < len(content):
        line = content[i]

        if line.startswith('def'):
            func_name = line.split()[1].split('(')[0]
            #func_name = line.strip('def').strip().strip(':')
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
                i, doc = parse_multi_line_doc(i, content)
            else:
                i, doc = parse_one_line_doc(i, content)

            func_dict[func_name] = doc
            in_func_scope = False

        elif in_func_scope and tri_quote not in line:
            # no docstring
            func_dict[func_name] = None
            in_func_scope = False
            i += 1

        if not in_func_scope:
            i += 1
            
    return func_dict


###
### start to parse .py files
###

module_dict = {}    # key: .py filename, value: dictionary with key is function name and value is function docstring

pyfilenames = [f for f in os.listdir('./') if f.endswith('py') and not f.startswith('_')]
for pyfname in pyfilenames:
    with open(pyfname, 'r') as pyf_content:
        content = pyf_content.readlines()
        func_dict = parse_pyfile(content)
        module_dict[pyfname] = func_dict
        

###
### write parsed result into README.md
###  
    
with open('README.md', 'w') as mdfile:
    write_content = []
    for module_name, sub_dict in module_dict.items():
        write_content.append(module_name)
        write_content.append('------')
        write_content.append('| Function | Description |')
        write_content.append('| -------- | ----------- |')
        
        for func_name, func_doc in sub_dict.items():
            #write_content.append(' - ' + func_name)
            #write_content.append('>' + func_doc)
            #write_content.append(f'| {func_name} | {func_doc} |')
            write_content.append(f'| <font color="#a77864"> **{func_name}** </font> | {func_doc} |')
        
        write_content.append('\n')
        write_content.append('******')
        
    write_string = '\n'.join(write_content)
    mdfile.write(write_string)