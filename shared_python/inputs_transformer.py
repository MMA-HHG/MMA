# Add standardised transformed inputs into hdf5-archive
import MMA_administration as MMA
import mynumerics as mn

def add_variables2hdf5(f,
                    global_input_names_to_jupyter_variables,
                    CUPRAD_names_to_jupyter_variables,
                    CTDSE_names_to_jupyter_variables,
                    CTDSE_outputs_to_jupyter_names,
                    list_of_CTDSE_outputs,
                    Hankel_names_to_jupyter_variables):

    # global inputs
    if not(global_input_names_to_jupyter_variables == None):
        global_inputs = f.create_group(MMA.paths['global_inputs'])
        for dset_name, (value, unit) in global_input_names_to_jupyter_variables.items():
            mn.adddataset(global_inputs,dset_name,value,unit)
    
    # CUPRAD
    if not(CUPRAD_names_to_jupyter_variables == None):
        CUPRAD_inps = f.create_group(MMA.paths['CUPRAD_inputs'])
        for dset_name, (value, unit) in CUPRAD_names_to_jupyter_variables.items():
            mn.adddataset(CUPRAD_inps,dset_name,value,unit)


    # CTDSE
    if not(CTDSE_names_to_jupyter_variables == None):
        # regular inputs
        CTDSE_inps = f.create_group(MMA.paths['CTDSE_inputs'])
        for dset_name, (value, unit) in CTDSE_names_to_jupyter_variables.items():
            mn.adddataset(CTDSE_inps,dset_name,value,unit)
    
        # CTDSE prints
        mn.adddataset(CTDSE_inps,'print_GS',1.,'[-]')                        
        for CTDSE_name, jupyter_name in CTDSE_outputs_to_jupyter_names.items():
            mn.adddataset(CTDSE_inps,CTDSE_name,int(jupyter_name in list_of_CTDSE_outputs),'[-]')

        # obsolete inputs: to be removed in a future release
        for dset_name in ('print_F_Efield_M2','print_F_Source_Term_M2','InterpByDTorNT'):
            mn.adddataset(CTDSE_inps,dset_name,0,'[-]')
        mn.adddataset(CTDSE_inps,'Ntinterp',1,'[-]')


    # Hankel
    if not(Hankel_names_to_jupyter_variables == None):
        Hankel_inps = f.create_group(MMA.paths['Hankel_inputs'])
        for dset_name, (value, unit) in Hankel_names_to_jupyter_variables.items():
            mn.adddataset(Hankel_inps,dset_name,value,unit)


def line_creator(name_dict, type_dict):
    content = ''
    for dset_name, (value, unit) in name_dict.items():
        if dset_name in type_dict['I']:
            content += dset_name + '\t' + str(value) + '\t' + 'I' + '\t' + unit + '\n'
        elif dset_name in type_dict['R']:
            content += dset_name + '\t' + str(value) + '\t' + 'R' + '\t' + unit + '\n'
        elif dset_name in type_dict['S']:
            content += dset_name + '\t' + value.decode() + '\t' + 'S' + '\t' + unit + '\n'
        elif dset_name in type_dict['R-array']:
            content += '$array' + '\t'  + dset_name + '\t' + 'R' + '\t' + unit + '\t' +\
                       ' '.join([str(foo) for foo in value]) +'\n'
        else:
            raise KeyError('The variable is not in the respective list (in type_dict).')

    return content

# transform to text formatted for universal inputs
def variables2text(global_input_names_to_jupyter_variables,
                   CUPRAD_names_to_jupyter_variables,
                   CTDSE_names_to_jupyter_variables,
                   CTDSE_outputs_to_jupyter_names,
                   list_of_CTDSE_outputs,
                   Hankel_names_to_jupyter_variables):

    # global inputs
    if not(global_input_names_to_jupyter_variables == None):
        content = '$change_group'+'\t'+ MMA.paths['global_inputs'] +'\n\n'
        content += line_creator(global_input_names_to_jupyter_variables, MMA.global_variable_type_lists)

    # CUPRAD inputs
    if not(CUPRAD_names_to_jupyter_variables == None):
        content += '\n\n' + '$change_group'+'\t'+ MMA.paths['CUPRAD_inputs'] +'\n\n'
        content += line_creator(CUPRAD_names_to_jupyter_variables, MMA.CUPRAD_variable_type_lists)

    # CTDSE inputs
    if not(CTDSE_names_to_jupyter_variables == None):
        content += '\n\n' + '$change_group'+'\t'+ MMA.paths['CTDSE_inputs'] +'\n\n'
        content += line_creator(CTDSE_names_to_jupyter_variables, MMA.CTDSE_variable_type_lists)
    
        # CTDSE prints
        content += 'print_GS' + '\t' + '1' + '\t' + 'I' + '\t' + '[-]' + '\n'
        for CTDSE_name, jupyter_name in CTDSE_outputs_to_jupyter_names.items():
            content += CTDSE_name + '\t' + str(int(jupyter_name in list_of_CTDSE_outputs)) + '\t' +\
                       'I' + '\t' + '[-]' + '\n'
            # mn.adddataset(CTDSE_inps,CTDSE_name,int(jupyter_name in list_of_CTDSE_outputs),'[-]')


    # Hankel inputs
    if not(Hankel_names_to_jupyter_variables == None):
        content += '\n\n' + '$change_group'+'\t'+ MMA.paths['Hankel_inputs'] +'\n\n'
        content += line_creator(Hankel_names_to_jupyter_variables, MMA.Hankel_variable_type_lists)

    # for dset_name, (value, unit) in global_input_names_to_jupyter_variables.items():
    #     content += dset_name + '\t' + str(value) + '\t' + unit + '\n'

    # content += '\n\n' + '$change_group'+'\t'+ MMA.paths['CUPRAD_inputs'] +'\n\n'
    # for dset_name, (value, unit) in CUPRAD_names_to_jupyter_variables.items():
    #      content += dset_name + '\t' + str(value) + '\t' + unit + '\n'

    return content