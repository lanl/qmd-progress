import re

import copy

def get_public(fortran_code):

   # Define the regex pattern to match lines starting with 'public ::' 
   pattern = re.compile(r'public ::(.*)')

   # Find matches in the fortran code
   matches = pattern.findall("".join(fortran_code))

   public_subroutines = []

   for match in matches:
      # Split by comma to get individual subroutines, strip whitespace
      subroutines = [subroutine.strip() for subroutine in match.split(',')]
      public_subroutines.extend(subroutines)

   # print the list of public subroutines
   print("Public subroutines are \n")
   for subroutine in public_subroutines:
       print(f"  public :: {subroutine}_c")
   print("")
   return public_subroutines


def transform_fortran_file(filename):

    #subroutine_re = re.compile(r'subroutine\s+([a-zA-Z0-9_]+)\s*\((.*)\)')
    subroutine_re = re.compile(r'subroutine\s+([a-zA-Z0-9_]+)\s*\(([\s\S]*)\)', re.DOTALL)

    #argument_re = re.compile(r'(real\(8\)|real\(dp\)|integer|type\(bml_matrix_t\))(, .+?)?\s+::\s+([a-zA-Z0-9_,\s]+)')
    #argument_re = re.compile(r'(real\(8\)|real\(dp\)||integer|type\(bml_matrix_t\))(, .+?)?\s+::\s+([a-zA-Z0-9_,\s\(\):]+)')
    argument_re = re.compile(r'(real\(8\)|real\(dp\)|integer|logical|type\(bml_matrix_t\)|character\(\d+\)|character\(len=\*\))(, .+?)?\s+::\s+([a-zA-Z0-9_,\s\(\):]+)')
    argument_re = re.compile(r'(real\(8\)|real\(dp\)|real\(PREC\)|integer|logical|type\(bml_matrix_t\)|character\(\d+\)|character\(len=\*\))(,.+?)?\s+::\s+([a-zA-Z0-9_,\s\(\):]+)')


    # Replacement rules for argument types
    type_replacement = {
        "real(8)": "real(c_double)",
        "real(dp)": "real(c_double)",
        "real(PREC)": "real(c_double)",
        "type(bml_matrix_t)": "type(c_ptr)",
        "integer": "integer(c_int)",
        "character(len=*)":"character(c_char)",
        "character(10)":"character(c_char)",
        "character(20)":"character(c_char)",
        "character(50)":"character(c_char)",
        "character(3)":"character(c_char)",
        "character(2)":"character(c_char)",
        "logical":"logical(c_bool)"
    }

    # Replacement rules for argument types (for c header)
    type_replacement_header = {
        "real(8)": "double",
        "real(dp)": "double",
        "real(PREC)": "double",
        "type(bml_matrix_t)": "bml_matrix_t*",
        "integer": "int",
        "character(len=*)":"char*",
        "character(10)":"char*",
        "character(20)":"char*",
        "character(50)":"char*",
        "character(3)":"char*",
        "character(2)":"char*",
        "logical":"int*"
    }

    with open(filename, "r") as f:
        content = f.readlines()

    public_subroutines = get_public(content)

    new_content = []
    header_content = []
    is_inside_subroutine = False
    argument_names = []
    argument_declarations = []
    argument_assignments = []
    subroutine_name = ""
    bml_matrix_t_args = []

    for k, line in enumerate(content):
        if "&" in line:
            line = line.replace("\n","").strip()
            line = line + content[k+1].strip()
            line = line.replace("&", "")
        match = subroutine_re.match(line.strip())
        if match:
            subroutine_name = match.group(1)
            argument_names = [arg.strip() for arg in match.group(2).split(',')]
            for iarg, argname in enumerate(argument_names):
                if "&" in argname:
                    print(f"debug-zy: {argname} has character &")
                    argument_names[iarg] = argname.split("&")[0]


            ispublic = subroutine_name in public_subroutines
            print(f"\ndebug-zy: a subroutine {subroutine_name} is found!", ispublic)
            print(f"debug-zy: its arguments are: {argument_names}")
            # skip the subroutine if it's not public
            if not ispublic: continue
            is_inside_subroutine = True
            arg2typ = {}
            continue

        if is_inside_subroutine:
            if 'end subroutine' in line:
                argument_names_c = [arg+"_c" if arg in bml_matrix_t_args else arg for arg in argument_names]
                new_subroutine_header = f'  subroutine {subroutine_name}_c({", ".join(argument_names_c)}) bind(C, name="{subroutine_name}")'
                new_subroutine_body = '\n'.join(argument_declarations + argument_assignments)
                new_subroutine_call = f'    call {subroutine_name}({", ".join(argument_names)})'
                new_content.append(new_subroutine_header)
                new_content.append(new_subroutine_body)
                new_content.append(new_subroutine_call)

                new_content.append('  end subroutine ' + subroutine_name + '_c\n')
                tmpheader = []
                for arg in argument_names:
                    if arg in arg2typ:
                        #print(f"{arg} type is {arg2typ[arg]}, {type_replacement_header[arg2typ[arg]]}")
                        tmpheader.append(f'{type_replacement_header[arg2typ[arg]]} {arg}')
                header_content.append(f'\nvoid {subroutine_name}({", ".join(tmpheader)});')
                #header_content.append(f'void {subroutine_name}({", ".join(type_replacement.get(arg, arg)+" "+arg+"_c" if arg in bml_matrix_t_args else arg for arg in argument_names_c)});')
                #header_content.append(f'void {subroutine_name}({", ".join(type_replacement[arg_type]+" "+arg for arg_type, arg in argument_names_c)});')
                #header_content.append( f'void {subroutine_name}({", ".join(argument_names_c)});')

                is_inside_subroutine = False
                argument_names = []
                argument_declarations = []
                argument_assignments = []
                bml_matrix_t_args = []
                continue

            match = argument_re.match(line.strip())
            if match:
                arg_type = match.group(1)
                arg_modifiers = match.group(2)
                print("matched group(3)=", match.group(3))
                #arg_names = [name.strip() for name in match.group(3).split(",")]

                # Remove any spaces within array declarations
                input_string = match.group(3)
                # Split the input string by comma using lookahead assertion
                arg_names = re.split(r',(?![^(]*\))', input_string)
                arg_names = [arg.strip() for arg in arg_names]
                print("debug: arg_names=", arg_names)

                """
                # remove (:), (:,:) from the arguments
                for i in range(len(arg_names)):
                    # Remove (:,:) and (:)
                    if "(:,:)" in arg_names[i]:
                        arg_names[i] = arg_names[i].replace("(:,:)", "").strip()
                    if "(:)" in arg_names[i]:
                        arg_names[i] = arg_names[i].replace('(:)', '').strip()
                    if "(:" in arg_names[i]:
                        arg_names[i] = arg_names[i].replace("(:", "").strip()
                """

                # create arg to type mapping for c header
                for arg_name in arg_names:
                    if "(" in arg_name:
                        arg_name = arg_name.split("(")[0]
                    print(f"debug-zy: arg_name={arg_name}, type is {arg_type}")
                    if arg_name in argument_names:
                        #print(f"debug-zy: arg_name={arg_name}, type is {arg_type}")
                        arg2typ[arg_name] = arg_type

                if arg_type in type_replacement:
                    new_arg_type = type_replacement[arg_type]
                    print("new arg_type is:", new_arg_type)
                    if arg_modifiers is None:
                        arg_modifiers = ", value"
                    else:
                        arg_modifiers = arg_modifiers.replace("intent(inout)", "value")
                        arg_modifiers = arg_modifiers.replace("intent(in)", "value")
                        if "allocatable" in arg_modifiers or "optional" in arg_modifiers:
                            arg_modifiers = arg_modifiers.replace("value", "")

                    for arg_name in arg_names:
                        tmparg_modifiers = copy.copy(arg_modifiers)

                        #index = arg_name.find('(')
                        ## Extract the substring before the opening parenthesis
                        #tmparg_name = arg_name[:index]

                        if "(:,:)" in arg_name:
                            tmparg_name = arg_name[:-5]
                            tmparg_modifiers = tmparg_modifiers.replace("value", "target")
                        elif "(:)" in arg_name:
                            tmparg_name = arg_name[:-3]
                            tmparg_modifiers = tmparg_modifiers.replace("value", "target")
                        elif "(" in arg_name:
                            tmparg_name = arg_name.split("(")[0]
                            tmparg_modifiers = tmparg_modifiers.replace("value", "target")
                        else:
                            tmparg_name = arg_name
                        print("debug: arg_name/tmparg_name=", arg_name, tmparg_name)
                        print("debug: argument_names=", argument_names)

                        if tmparg_name in argument_names:
                        #if arg_name in argument_names:
                            print(f"{tmparg_name} is added to the declaration\n")
                            if arg_type == "type(bml_matrix_t)":
                                bml_matrix_t_args.append(arg_name)
                                argument_declarations.append(f'    {new_arg_type}{tmparg_modifiers} :: {arg_name}_c')
                                argument_declarations.append(f'    {arg_type} :: {arg_name}')
                                argument_assignments.append(f'    {arg_name}%ptr = {arg_name}_c')
                            else:
                                argument_declarations.append(f'    {new_arg_type}{tmparg_modifiers} :: {arg_name}')
                        else:
                            print(f"{tmparg_name} is not added to the declaration\n")

    new_filename = filename[:-4] + "_c.F90"
    print("c_wrappers are:")
    with open(new_filename, "w") as f:
        print('\n'.join(new_content))
        f.write('\n'.join(new_content))

    print("c headers are:")
    #header_filename = filename[:-4] + ".h"
    header_filename = "prg_progress_mod.h"
    with open(header_filename, "a") as f:
        print('\n'.join(header_content))
        f.write('\n'.join(header_content))
    f.close()

# Example usage
import glob, sys

files=glob.glob(f"src/*mod.F90")
for fname in files:
    print(fname)
    #fname = "prg_densitymatrix_mod.F90"
    transform_fortran_file(fname)
    #if "prg_ewald_mod" in fname: sys.exit()
    #if "prg_graphsolver_mod" in fname: sys.exit()
