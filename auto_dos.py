#!/usr/bin/env python
##########################################################
# auto_dos_new.py Script 0.3
#
# coded by John McLeod
# 2011 03 01
#
# This script is designed to generate and sum the DOS from all atoms of a given
# type in a structure
#
# x lawp -qtl MUST have already been run for this script to work
#
# the directory must be a "proper" w2web directory, i.e. the name of the folder should
# be the same as the <case> for the calculation
#
# this script requires a command line argument specifying the directory for the calculation
#
# optionally, this script may be called with the additional arguments:
# 	-All	specifies a distinct DOS file to be created for each atom, 
#		rather than summing by atomic species
#	-NaN	specifies that NaN entries should be recalculated as linear
#		values between the available endpoints
#	-Sym	specifies that all the symmetries should be created, rather than just s, p, and d
#	-Tot	specifies that the total DOS should also be calculated
#	-X#	specifies that the unbroadened XES should be calculated
#		if # is K, then only K-edges will be calculated
#		if # is L, then only L2,3-edges will be calculated
#		if # is M, then only M4,5-edges will be calculated
#		if # is not specified, then all available edges will be calculated
#	-brd	specifieds that the XES should also be broadened.
#		this uses a fixed spectrometer resolution of 1e3
#		and lorentz broadening with a scaling factor of 3 eV
#		and a core-hole lifetime parameter using Auger yield for isolated atoms
#		(see Schwarz, et. al., J. Phys. F: Metal Phys., 9(12), 2509 (1979).)
##########################################################
# CHANGELOG
# 0.32
# fix to handle multiplicity higher than 9
# 0.31
# fixed crash for non -up/-dn cases
#
# 0.3
# implemented -up and -dn
# fixed NPT= being appended to filenames
# support for more than 99 atoms in a structure
# fixed deprecation warning in range
##########################################################


# top level imports
import sys, os, subprocess, math
from scipy.interpolate import interp1d


##########################################################
# function get_params
# 
# this is a generic function designed to open a particular file,
# identify lines containing a particular FLAG, and retrieve
# a set of values from that line
#
# optionally, a FLAG_IDX can be specified identifying the location of FLAG on the line
#		(if FLAG_IDX is given, then FLAG can be a substring of that item)
#		(if FLAG_IDX is NOT given, then FLAG must be the entire string)
# optionally, a DELIMITER can be specified identifying the character that separates data values
# 		(if not specified, blank space is assumed)
# optionally, a IDX_LEN can be specified for identifying the number of characters expected for the IDX value
#		(if not specified, no length is assumed)
#		(this is used for string-related values that might have a space in them)
##########################################################
def get_params(file_name, flag, idx_list, flag_idx=-1, delimiter='', idx_len=[]):
    # first make sure that this file exists
    if os.path.exists(file_name) == False:
        print("ERROR: ", file_name, "does not exist! Did you forget -up or -dn?")
        return -1

    # initialize an empty tuple for the return data
    return_data = ()

    # open the file
    file_data = open(file_name, 'r')
    # scan through the file
    for line in file_data:
        # make sure this isn't a blank line
        if (line == '') or (line == '\n'):
            continue

        # if ( flag_idx == 0) :
        #	print line[:5]+' '+line[6:]
        # split up the line
        if delimiter == '':
            data = line.split()
        else:
            data = list(my_split(line, delimiter))

        # check to see if this is a valid line
        # line is valid if the <FLAG_IDX>th value of DATA is FLAG,
        # or if an element of DATA is FLAG (if no FLAG_IDX is specified)
        if ((flag_idx != -1) and (data[flag_idx].find(flag) != -1)) or ((flag_idx == -1) and (data.count(flag) == 1)):
            # check if flag was a substring of any split element and re-split
            # (assumes flag is first occurrence in the string
            if (len(data[data[flag_idx].find(flag)]) > len(flag)):
                templine = line[:line.find(flag)] + ' ' + flag + ' ' + line[line.find(flag) + len(flag):]
                if delimiter == '':
                    data = templine.split()
                else:
                    data = my_split(templine, delimiter)

            # add data values to the tuple
            if len(return_data) == 0:
                for idx in idx_list:
                    # check to see if the length of the data value agrees with the specified length
                    if (idx_len == []) or (len(idx_len) < len(idx_list)) or (len(data[idx]) == idx_len[idx]):
                        return_data += [data[idx]],
                    # if it is too long, only return the first few characters
                    elif len(data[idx]) > idx_len[idx]:
                        return_data += [data[idx][:idx_len[idx]]],
                    # if it is too short, merge the two entries
                    # (this will probably be the most common operation if value length is imporant, since spaces
                    #	in string values may cause that string to be separated into several different data elements
                    elif len(data) > idx + 1:
                        return_data += [(data[idx] + data[idx + 1])],
            else:
                for i, idx in enumerate(idx_list):
                    # check to see if the length of the data value agrees with the specified length
                    if (idx_len == []) or (len(idx_len) < len(idx_list)) or (len(data[idx]) == idx_len[idx]):
                        return_data[i].append(data[idx])
                    # if it is too long, only return the first few characters
                    elif len(data[idx]) > idx_len[idx]:
                        # return_data[i].append( data[idx][:idx_len[idx]] )
                        return_data[i].append(data[idx])
                    # if it is too short, merge the two entries
                    # (this will probably be the most common operation if value length is imporant, since spaces
                    #	in string values may cause that string to be separated into several different data elements
                    elif len(data) > idx + 1:
                        return_data[i].append(data[idx] + data[idx + 1])

    file_data.close()

    return return_data


join = lambda x: sum(x, [])


# function my_split
# iterate through seps to split string s
def my_split(s, seps):
    fragments = [s]
    for token in seps:
        fragments = join(f.split(token) for f in fragments)
    fragments = filter(None, fragments)
    return fragments


##########################################################
# function make_dos
# 
# this is a function designed to write a TETRA input file
# and use the shell command 'x tetra' to generate DOS files
#
# this function requires a CASE string for the WIEN2k <case> name
# this function requires an ATOM string for the WIEN2k atom number
# this function requires a SYMMETRY string of comma-separated symmetry names
# this functional optionally requires a FLAG to specify whether full symmetry (flag = 1)
#	or only s, p, and d symmetry (flag = -1) are generated
##########################################################
def make_dos(case, atom, symmetry, flag=-1):
    # open a TETRA input file
    int_file = open(case + '.int', 'w')
    # write header information
    int_file.write("DOS parameters generated by auto_dos.py\n")
    int_file.write(" -1.50 0.002 1.50 0.003\n")

    # determine how many symmetries we need for this atom
    if (flag == -1) and (len(symmetry.split(',')) > 3):
        # just use s, p, d
        dos_num = '3'
    else:
        # use all
        dos_num = str(len(symmetry.split(',')))
        flag = 1
    # write the number of DOS cases and finish the header information
    int_file.write("    " + dos_num + "    N   0.000\n")

    # write the symmetries into the input file
    for idx, sym in enumerate(symmetry.split(',')):
        temp_flag = 1
        if sym == '0':
            sym = 's'
            temp_flag = flag
        elif sym == '1':
            sym = 'p'
            temp_flag = flag
        elif sym == '2':
            sym = 'd'
            temp_flag = flag
        if temp_flag == flag:
            int_file.write('   ' + atom + '   ' + str(idx + 1) + '    ' + sym + '\n')
    int_file.close()

    # run TETRA to generate the DOS
    # (with hacked conditional operator to add -up or -dn)

    if ((updn == "up") | (updn == "dn")):
        os_call = subprocess.call(['x', 'tetra', '-' + updn])
    else:
        os_call = subprocess.call(['x', 'tetra'])

    if os_call != 0:
        print("Running TETRA failed!")

    # BUG TESTING ##################
    # copy <case>.int for checking #
    # os_call = subprocess.call( ['cp', case + '.int', case + '.int_' + atom] )
    ################################

    # return the TETRA run code
    return os_call


##########################################################
# function make_xes
# 
# this is a function designed to write a XSPEC input file
# and use the shell command 'x xspec' to generate XES files
#
# this function requires a CASE string for the WIEN2k <case> name
# this function requires an ATOM string for the WIEN2k atom number
# this functional optionally requires a FLAG to specify whether a particular
#	XES symmetry (1, 2, 3 for K, L2,3, and M4,5) are generated
#	default is the K-XES
##########################################################
def make_xes(case, atom, flag=1):
    # open a TETRA input file
    int_file = open(case + '.inxs', 'w')
    # write header information
    int_file.write("XES parameters generated by auto_dos.py\n")
    # write the atom number
    int_file.write(' ' + atom + "          (atom)\n")

    # make sure we have a valid flag
    if (int(flag) < 1) or (int(flag) > 3):
        print("Must be a valid edge (K=1, L2,3=2 or M4,5=3), unknown edge specified : ", flag)
        return 1

    # write the edge information
    int_file.write(str(flag) + "          (n core)\n")
    int_file.write(' ' + str(flag - 1) + "          (l core)\n")

    # write generic information
    int_file.write(' 0, 0.5, 0.5\n')
    int_file.write(' -25, 0.1, 5\n')
    int_file.write('EMIS\n')
    # since we broaden on our own, these parameters are not essential
    int_file.write('0.3           (S)\n')
    int_file.write('0.3           (gamma0)\n')
    int_file.write('0.8           (W)\n')
    int_file.write('AUTO\n')
    int_file.close()

    # run XSPEC to generate the DOS

    if ((updn == "up") | (updn == "dn")):
        os_call = subprocess.call(['x', 'xspec', '-' + updn])
    else:
        os_call = subprocess.call(['x', 'xspec'])

    if os_call != 0:
        print("Running XSPEC failed!")
        return os_call

    # move the XSPEC file in case we are generating more
    os_call = subprocess.call(['mv', case + '.txspec' + updn, case + '.txspec' + updn + '_' + 'KLM'[flag - 1]])
    if os_call != 0:
        print("Couldn't move .txspec file!")

    # read through the created dos1ev file to find the band edges
    dos_file = open(case + '.dos1ev' + updn, 'r')
    # initialize energy and dos arrays
    e = []
    dos = []
    # cycle through the file
    for line in dos_file:
        # discard all header lines
        if line[0] == '#':
            continue
        # ignore all DOS at energies greater than the fermi level
        if float(line.split()[0]) > 0.0:
            continue
        # add the values to the list
        e.append(float(line.split()[0]))
        dos.append(float(line.split()[-1]))
    dos_file.close()

    # initialize array for storing the energies of the band edges
    e_edge = []
    band_flag = 0
    # cycle through the energies
    for idx in range(len(e_edge)):
        # look at the dos starting at the fermi edge and going backwards
        # if the DOS is greater than zero, this is an energy band
        if (dos[-idx - 1] > 0):
            # if we haven't picked the top edge of this band yet, grab it
            if (band_flag == 0):
                e_edge.append(e[-idx - 1])
                band_flag = 1
        else:
            # if we haven't picked the bottom edge of this band yet, grab it
            if (band_flag == 1):
                e_edge.append(e[-idx - 1])
                band_flag = 0
    # make sure we have a lower edge to this band
    if len(e_edge) == 1:
        e_edge.append(e[0])

    return e_edge


##########################################################
# function make_xas
# 
# this is a function designed to write a XSPEC input file
# and use the shell command 'x xspec' to generate XAS files
#
# this function requires a CASE string for the WIEN2k <case> name
# this function requires an ATOM string for the WIEN2k atom number
# this functional optionally requires a FLAG to specify whether a particular
#	XAS symmetry (1, 2, 3 for K, L2,3, and M4,5) are generated
#	default is the K-XAS
##########################################################
def make_xas(case, atom, flag=1):
    # open a TETRA input file
    int_file = open(case + '.inxs', 'w')
    # write header information
    int_file.write("XES parameters generated by auto_dos.py\n")
    # write the atom number
    int_file.write(' ' + atom + "          (atom)\n")

    # make sure we have a valid flag
    if (int(flag) < 1) or (int(flag) > 3):
        print("Must be a valid edge (K=1, L2,3=2 or M4,5=3), unknown edge specified : ", flag)
        return 1

    # write the edge information
    int_file.write(str(flag) + "          (n core)\n")
    int_file.write(' ' + str(flag - 1) + "          (l core)\n")

    # write generic information
    int_file.write(' 0, 0.5, 0.5\n')
    int_file.write(' -5, 0.1, 45\n')
    int_file.write('ABS\n')
    # since we broaden on our own, these parameters are not essential
    int_file.write('0.3           (S)\n')
    int_file.write('0.3           (gamma0)\n')
    int_file.write('0.8           (W)\n')
    int_file.write('AUTO\n')
    int_file.close()

    # run XSPEC to generate the DOS
    if ((updn == "up") | (updn == "dn")):
        os_call = subprocess.call(['x', 'xspec', '-' + updn])
    else:
        os_call = subprocess.call(['x', 'xspec'])
    if os_call != 0:
        print("Running XSPEC failed!")
        return os_call

    # move the XSPEC file in case we are generating more
    os_call = subprocess.call(['mv', case + '.txspec' + updn, case + '.txspec' + updn + '_' + 'KLM'[flag - 1]])
    if os_call != 0:
        print("Couldn't move .txspec file!")

    # return the edges of the XAS, no fancy band edge stuff here
    return [0.0, 25.0]


##########################################################
# function clean_dos
# 
# this is a function designed to clean up DOS?EV files
#
# this function requires a CASE string for the WIEN2k <case> name
# this function optionally requires an ATOM string for the new file name
# this function optionally requires a NAN flag to identify whether or not to discard NaN values
#	(defaults to including NaN values)
# this function merges all DOS?EV files, writes clean headers, and separates data by commas
#
# this function produces a CASE.DOS_ATOM file or a CASE.DOS_new file if
# 	ATOM is not specified
##########################################################
def clean_dos(case, atom='new', NaN=-1):
    # check to see if DOS file exists
    if os.path.exists(case_name + '.dos1ev' + updn) == False:
        print("No DOS1EV file found!")
        return -1

    data = []
    energy = []
    symmetry = ''

    # read in the dos?ev files if they exist
    for stem in ['.dos1ev' + updn, '.dos2ev' + updn, '.dos3ev' + updn]:
        if os.path.exists(case_name + stem) == False:
            continue

        # open the file
        dos_data = open(case_name + stem, 'r')

        # cycle through the data
        for line in dos_data:
            # get symmetry headers
            if line[0] == '#':
                # symmetry headers are on line starting with '# ENERGY'
                if line[:8] == '# ENERGY':
                    # add everything except initial hash and final '\n' character
                    if symmetry == '':
                        symmetry += line[1:-1]
                    else:
                        # don't add the 'ENERGY' tag twice
                        symmetry += line[1:-1].split('ENERGY')[-1]
                continue

            # get the energy of these DOS data
            e = line.split()[0]

            # check to see if this is the first set of data we've encountered
            # if so, add it, removing '\n' character at end
            #	(even if '\n' character isn't there we will just remove the least significant decimal, no big deal)
            if data == []:
                data.append(line[:-1])
                energy.append(e)
            # if this energy value is greater than the largest in the previous data, add this DOS data to the end
            elif float(e) > float(energy[-1]):
                data.append(line[:-1])
                energy.append(e)
            # otherwise, append DOS data to previous entry
            elif energy.count(e) == 1:
                entries = line.split()
                for item in entries[1:]:
                    data[energy.index(e)] += '   ' + item
            # if energy doesn't exist in the list, we may have a problem
            else:
                print("Energy ", e, " in file ", stem, " does not exist in list!")

        dos_data.close()

    # make sure we have some data
    if data == []:
        print("No DOS data acquired!")
        return -1

    # check to see if the cleaned DOS file already exists
    if os.path.exists(case_name + '.DOS_' + atom):
        print("WARNING! We are overwriting clean DOS file: ", case_name + '.DOS_' + atom)

    # write a new dos file
    dos_data = open(case_name + '.DOS_' + atom, 'w')
    # write the column headers
    for val in symmetry.split()[:-1]:
        if (val != 'ENERGY'):
            dos_data.write(val[-1]+', ') #strip off the atom number
        else:
            dos_data.write(val + ', ')
    dos_data.write(symmetry.split()[-1] + '\n')

    # write the DOS data as comma-separated values
    for item in data:
        # ignore data if it contains a NaN value and the NAN flag is 1
        if (NaN == 1) and (item.find('NaN') != '-1'):
            continue
        # otherwise add the data to the file
        for val in item.split()[:-1]:
            dos_data.write(val + ', ')
        dos_data.write(item.split()[-1] + '\n')

    dos_data.close()

    # return the file name
    return case_name + '.DOS_' + atom


##########################################################
# function read_txspec
# 
# this is a function designed to read a TXSPEC file
#
# this function requires a file name
#
# this function returns the x and total unbroadened xspec data
##########################################################
def read_txspec(file_name=''):
    # check to make sure this file exists
    if os.path.exists(file_name) == False:
        print("TXSPEC file ", file_name, " was not found!")
        return -1

    # open the file
    tx_file = open(file_name, 'r')
    # initialize lists
    e = []
    xes = []

    # read through the data
    for line in tx_file:
        # this should always be a line of actual data
        if line.strip(' -+.,0123456789eE\n') != '':
            print("Are you sure this is a TXSPEC file?")
            continue
        # add the values to the list
        if (line != '') and (line != '\n'):
            e.append(float(line.split()[0]))
            xes.append(float(line.split()[1]))

    # make sure we have enough data
    # in the case of XAS, the calculation can sometimes run out of states
    # before the end of the file
    if len(e) < 301:
        for idx in range(301 - len(e)):
            e.append(e[-1] + 0.1)
            xes.append(xes[-1])

    return e, xes


##########################################################
# function broaden_txspec
# 
# this is a function designed to broaden TXSPEC data
#
# this function requires a list of energies and unbroadened transition
#	probabilities
#
# optionally this function requires a list of band edge energies
#	to apply variable broadening
#	(default is to broaden based on entire energy range)
# optionally, this function requires a life-time broadening energy
#	(if not provided a default of 0.2 eV is used)
# optionally, this function requires a spectrometer resolution
#	(if not provided a default of 0.1 eV is used)
# optionally this function requires an energy splitting between edges
#	(such as L2,3)
#	(default is to apply no splitting)
# optionally, this function requires a ratio between split edges
#	(default is to scale edges equally if an energy splitting is provided)
#	(if no energy splitting is provided this value is ignored)
# optionally, this function requires a W scale factor
#	(if not provided a default of 0.8 eV is used)
# optionally, this function requires an energy shift value
#	(if provided, then the spectrometer resolution should be the RESOLVING POWER, e.g. 1000)
#	(gaussian resolution will then be (shift + energy) / res)
#	(otherwise, a default of -1.0 is used and `res' is used as full broadening width
# this function returns the x and total unbroadened xspec data
##########################################################
def broaden_txspec(energy, xes, edges=[], gamma0=0.1, res=0.1, split=-1.0, ratio=1.0, W=0.4, shift=-1.0):
    # figure out if a list of edges were supplied, if not just use fermi level and end of spectrum
    if len(edges) < 2:
        edges = [0.0, energy[0]]

    # check to see that data lengths are consistent
    if len(energy) != len(xes):
        print("Incorrectly formatted Energies and XES data!")
        return -1

    # figure out energy spacing
    delta_e = energy[1] - energy[0]
    #	print ('##### ENERGY STEP #####', delta_e)
    #	print ('##### ENERGY SPLIT #####', split)

    # create empty broadened XES data
    xes_brd = [0.0 for x in xes]

    # cycle through data to apply life-time broadening
    for idx, e1 in enumerate(energy):
        for jdx, e2 in enumerate(energy):
            # find band edges :
            e_upper = 0.0
            e_lower = energy[0]
            n_upper = 0
            n_lower = 1
            for kdx, e_edge in enumerate(edges):
                # if the edge is higher than the current energy
                # we are below that upper band edge
                if e2 <= e_edge:
                    e_upper = e_edge
                    n_upper = kdx
                # if the edge is lower than the current energy
                # we are above that lower band edge
                elif e2 > e_edge:
                    e_lower = e_edge
                    n_lower = kdx

            # make sure we haven't royally fucked up
            if e_upper == e_lower:
                print("The top of the band is the same as the bottom of the band!")
                print(edges, e_upper, e_lower)

            # determine the gamma factor
            gamma = gamma0 / 2.0
            # are we inside a populated band? If so, n_upper will be even and n_lower will be odd
            if (n_upper % 2 == 0) and (n_lower % 2 == 1):
                gamma += W * (1.0 - (e2 - e_upper) / (e_upper - e_lower)) ** 2
                # add extra W factors for deep bands
                gamma += (n_upper / 2) * W
            # if we are not in a band I don't think there should be any states, but just in case
            else:
                gamma += W

            # apply lorentz broadening
            xes_brd[idx] += xes[jdx] / math.pi * (
                    math.atan((e2 - e1 + delta_e) / gamma) - math.atan((e2 - e1 - delta_e) / gamma))

    # split and scale the spectra, if appropriate
    if split != -1.0:
        # get the approximate split in indices of the spectra
        idx_split = int(split / delta_e)
        for idx in range(len(xes_brd) - idx_split):
            xes_brd[idx + idx_split] += ratio * xes_brd[idx]

    # apply spectrometer broadening, if appropriate
    if res > 0.0:
        # initialize a blank list
        xes_final = [0 for x in xes_brd]
        # cycle through data to apply gaussian instrumental broadening
        for idx, e1 in enumerate(energy):
            for jdx, e2 in enumerate(energy):
                # calculate the broadening width
                if shift > 0:
                    sigma = abs(float((e2 + shift) / res) / (2.0 * math.log(2.0)))
                else:
                    sigma = abs(res / (2.0 * math.log(2.0)))
                # apply gaussian broadening
                xes_final[idx] += xes_brd[jdx] * 2.0 * delta_e / math.sqrt(2.0 * math.pi) / sigma * math.exp(
                    -(e1 - e2) ** 2 / (2.0 * sigma ** 2))
    else:
        # otherwise just pass back the lorentz-broadened array
        xes_final = xes_brd

    # return the broadened XES
    return xes_final


##########################################################
# function clean_xes
# 
# this is a function designed to clean up TXSPEC files
#
# this function requires a CASE string for the WIEN2k <case> name
# this function optionally requires an ATOM string for the new file name
# this function optionally requires an SPEC string for the file type
#	default is XES
# this function optionally requires an E_BIND list of binding energies
#	to shift the XES/XAS files to the proper energy
# this function optionally requires an E_EDGE list of valence bands
#	to apply broadening
# this function optionally requires an ATOM_Z for the Z-number of the atom
#	this is also used for broadening
#
# this function produces a CASE.XES_ATOM file or a CASE.DOS_new file if
# 	ATOM is not specified
# all possible edges (K, L, M) will be included in this file if they are
# 	present
##########################################################
def clean_xes(case, atom='new', spec='XES', e_bind=[], e_edge=[], atom_Z=-1):
    # create blank tuple for data
    xes_data = ()
    header = ''
    tuple_len = 0

    if e_edge == None:
        return -1

    # determine resolution
    if spec == 'XAS':
        resolution = 0.2
    else:
        #		resolution = 0.6
        resolution = 0.4

    # check to see if a K-edge TXSPEC file exists
    if os.path.exists(case_name + '.txspec' + updn + '_K') == True:
        # get K-edge data
        (energy, txspec) = read_txspec(case_name + '.txspec' + updn + '_K')
        # if there is no band edges, just make the entire range the valence band
        if e_edge == []:
            e_edge = [energy[0], energy[-1]]
        # broaden K-edge data
        xspec = broaden_txspec(energy, txspec, e_edge, res=resolution)
        # shift the energy
        if e_bind == []:
            e_shift = energy
        else:
            e_shift = [e - e_bind[0] for e in energy]
        xes_data += (e_shift, txspec, xspec)
        tuple_len += 1
        # write header data
        header += 'K energy, K raw, K broadened'
    # check to see if a L-edge TXSPEC file exists
    if os.path.exists(case_name + '.txspec' + updn + '_L') == True:
        # get L-edge data
        (energy, txspec) = read_txspec(case_name + '.txspec' + updn + '_L')
        # if there is no band edges, just make the entire range the valence band
        if e_edge == []:
            e_edge = [energy[0], energy[-1]]
        # shift the energy and get splitting
        if len(e_bind) < 3:
            e_shift = energy
            e_split = -1
        else:
            e_shift = [e - e_bind[1] for e in energy]
            e_split = e_bind[2] - e_bind[1]

        # broaden K-edge data, account for L2,3 splitting and ideal ratio
        xspec = broaden_txspec(energy, txspec, e_edge, atom_Z, res=resolution, split=e_split, ratio=0.5)
        xes_data += (e_shift, txspec, xspec,)
        tuple_len += 1
        # write header data
        if header != '':
            header += ', '
        header += 'L23 energy, L23 raw, L23 broadened'
    # check to see if a M-edge TXSPEC file exists
    if os.path.exists(case_name + '.txspec' + updn + '_M') == True:
        # get M-edge data
        (energy, txspec) = read_txspec(case_name + '.txspec' + updn + '_M')
        # if there is no band edges, just make the entire range the valence band
        if e_edge == []:
            e_edge = [energy[0], energy[-1]]
        # shift the energy and get splitting
        if len(e_bind) < 5:
            e_shift = energy
            e_split = -1
        else:
            e_shift = [e - e_bind[3] for e in energy]
            e_split = e_bind[4] - e_bind[3]

        # broaden M-edge data, account for M4,5 splitting and ideal ratio
        xspec = broaden_txspec(energy, txspec, e_edge, res=resolution, split=e_split, ratio=0.66666667)
        xes_data += (e_shift, txspec, xspec,)
        tuple_len += 1

        # write header data
        if header != '':
            header += ', '
        header += 'M45 energy, M45 raw, M45 broadened'

    # check to see if a generic TXSPEC file exists
    if os.path.exists(case_name + '.txspec' + updn) == True:
        # get M-edge data
        (energy, txspec) = read_txspec(case_name + '.txspec' + updn)
        # if there is no band edges, just make the entire range the valence band
        if e_edge == []:
            e_edge = [energy[0], energy[-1]]
        # shift the energy and get splitting, just use the first binding energy value
        if len(e_bind) < 1:
            e_shift = energy
        else:
            e_shift = [e - e_bind[0] for e in energy]

        # broaden data
        xspec = broaden_txspec(energy, txspec, e_edge, res=resolution)
        xes_data += (e_shift, txspec, xspec,)
        tuple_len += 1

        # write header data
        if header != '':
            header += ', '
        header += 'XX energy, XX raw, XX broadened'

    # check to see that we actually got some data
    if xes_data == ():
        print("No " + spec + " data acquired!")
        return -1

    # check to see if the cleaned DOS file already exists
    if os.path.exists(case_name + '.' + spec + '_' + atom):
        print("WARNING! We are overwriting clean XES file: ", case_name + '.' + spec + '_' + atom)

    # write a new dos file
    xes_file = open(case_name + '.' + spec + '_' + atom, 'w')
    # write the column headers
    xes_file.write(header + '\n')

    # get lengths of data - they should be all the same but just in case...
    data_lengths = []
    for idx in range(tuple_len):
        data_lengths.append(len(xes_data[3 * idx]))

    # write the XES data as comma-separated values
    for idx in range(max(data_lengths)):
        # add the data to the file
        for jdx in range(tuple_len):
            # add data separators where appropriate
            if jdx > 0:
                xes_file.write(', ')
            # make sure we still have data in this XES range
            if idx < data_lengths[jdx]:
                for kdx in range(3):
                    xes_file.write(str(xes_data[3 * jdx + kdx][idx]))
                    if kdx != 2:
                        xes_file.write(', ')
            # otherwise write blanks
            else:
                xes_file.write(', , ')
        xes_file.write('\n')

    # close the file
    xes_file.close()

    # return the file name
    return case_name + '.' + spec + '_' + atom


##########################################################
# function sum_dos
# 
# this is a function designed to sum cleaned DOS_xxx files
# for all atoms of the same species
#
# ONLY the s, p, and d orbitals are summed, since symmetry projections
#	aren't universal for al atom species
#
# this function requires a FILE_LIST list of string for the clean DOS_xxx data
# this function optionally requires a SUM_FILE string specifying the output file name
#	if omitted, default name 'summed_dos' is used
# this function optionally requires a MULT_LIST list of integers specifying the multiplicity of each DOS site
#	if omitted, multiplicity defaults to 1 for each site
##########################################################
def sum_dos(file_list, sum_file='summed_dos', mult_list=[]):
    # populate multiplicity, if necessary
    if mult_list == []:
        mult_list = [1 for item in file_list]
    # get total number of atoms
    mult_tot = float(sum(mult_list))

    # initialize value keeping track of the previous number of s,p,d symmetries
    old_spd_idx = -1

    # create empty lists
    energy = []
    dos = []
    spd_idx = []
    header = []

    # read through all files
    for idx, dos_file in enumerate(file_list):
        # open the dos file
        dos_data = open(dos_file, 'r')

        # initialize indices for s, p, and d data
        spd_idx = []

        # read through the data
        for line in dos_data:
            # check if line is header
            # this is stupid because versions below 2.6 don't have the same nifty `translate' method
            # the `strip' method is less fool-proof, I think
            if (sys.version_info[1] > 5):
                temp_line = line.translate({ord(c): None for c in '-,. 0123456789Nae\n'})
            else:
                temp_line = line.strip('-,. 0123456789Nae\n')
            if temp_line != '':
                # look for indices of main symmetries: s, p, and d
                for i, s in enumerate(['s', 'p', 'd']):
                    # grab this header
                    header = line[:-1]

                    if header.split(', ').count(s) == 1:
                        spd_idx.append(header.split(', ').index(s))
                    else:
                        # if we can't find an s, p, and d orbital we are probably in trouble
                        print("No", s, " DOS found in file: ", dos_file)
                # skip the rest of this loop, since it deals with adding the data
                continue

            # check if we have found valid indices for s, p, and d orbitals
            # otherwise, kill this function since something has obviously gone wrong
            # the header for the orbitals should have come first in this file
            if spd_idx == []:
                print("We can't find any valid DOS in this file!")
                dos_data.close()
                return -1
            # otherwise, check to see if we found the same number of s,p,d orbitals in this file
            # as in the previous
            elif (old_spd_idx != -1) and (len(spd_idx) != old_spd_idx):
                print("WARNING! We seem to be dealing with a different number of s, p, and d orbitals!")
            # this is a crude check, but better than nothing

            # now get the energy and dos data from this line
            e = float(line.split(',')[0])
            d = [float(line.split(',')[val]) * mult_list[idx] / mult_tot for val in spd_idx]

            # if energy doesn't exist in the master list, add this data
            if energy.count(e) == 0:
                energy.append(e)
                dos.append(d)
            # otherwise add it to the existing data
            else:
                for i in range(len(d)):
                    dos[energy.index(e)][i] += d[i]

        # close the file
        dos_data.close()

        # save the spd_idx for consistency checking
        old_spd_idx = len(spd_idx)

    # open the summed file
    dos_data = open(sum_file, 'w')
    # write the header
    dos_data.write(header.split(', ')[0])
    for idx in spd_idx:
        dos_data.write(', ' + header.split(', ')[idx])
    dos_data.write('\n')
    # write in the data
    for idx in range(len(energy)):
        dos_data.write(str(energy[idx]))
        for val in dos[idx]:
            dos_data.write(', ' + str(val))
        dos_data.write('\n')

    dos_data.close()
    return 0


##########################################################
# function sum_xes
# 
# this is a function designed to sum cleaned XES_xxx files
# for all atoms of the same species
#
# this function requires a FILE_LIST list of string for the clean XES_xxx data
# this function optionally requires a SUM_FILE string specifying the output file name
#	if omitted, default name 'summed_dos' is used
# this function optionally requires a MULT_LIST list of integers specifying the multiplicity of each DOS site
#	if omitted, multiplicity defaults to 1 for each site
##########################################################
def sum_xes(file_list, sum_file='summed_xes', mult_list=[]):
    # populate multiplicity, if necessary
    if mult_list == []:
        mult_list = [1 for item in file_list]
    # get total number of atoms
    mult_tot = float(sum(mult_list))

    # create empty lists
    energy = []
    xes = []

    e_min = [0.0, 0.0, 0.0]
    e_max = [10000000.0, 10000000.0, 10000000.0]

    # get energy range quickly
    # read through all the files
    for xes_file in file_list:
        if (os.path.exists(xes_file) == False):
            return -1
        xes_data = open(xes_file, 'r')
        # cycle through all lines in the file
        for line in xes_data:
            # look for the first non-header line
            if line.strip('-,. 0123456789Nae\n') == '':
                for jdx in range(int(len(line.split(',')) / 3)):
                    if float(line.split(',')[3 * jdx]) > e_min[jdx]:
                        e_min[jdx] = float(line.split(',')[3 * jdx])
                    if float(line.split(',')[3 * jdx]) < e_max[jdx]:
                        e_max[jdx] = float(line.split(',')[3 * jdx])
                xes_data.close()
                break

    # read through all files
    for idx, xes_file in enumerate(file_list):
        # open the dos file
        xes_data = open(xes_file, 'r')

        # create empty lists for temporary data storage
        e = []
        x = []

        # read through the data
        for line in xes_data:
            # check if line is header
            # this is stupid because versions below 2.6 don't have the same nifty `translate' method
            # the `strip' method is less fool-proof, I think
            if (sys.version_info[1] > 5):
                temp_line = line.translate(None, '-,. 0123456789Nae\n')
            else:
                temp_line = line.strip('-,. 0123456789Nae\n')
            if temp_line != '':
                # skip the rest of this loop, since it deals with adding the data
                continue

            # now get the energy and dos data from this line
            line_data = line.split(',')
            # there may be more than one energy, since there may be more than one spectra
            # these energies are the even multiples of 3, though (starting with column 0)
            e.append([float(line_data[3 * i]) for i in range(int(len(line_data) / 3))])
            x.append([float(line_data[3 * i + 2]) for i in range(int(len(line_data) / 3))])

        # close the file
        xes_data.close()

        # check if master energy list exists, if not create it
        if energy == []:
            for jdx in range(len(e[0])):
                # energy spacing
                de = 0.1
                min_e = round(e_min[jdx] + 0.05, 1)
                max_e = round(e_max[jdx] + 29.95, 1)
                n_jdx = (max_e - min_e) / de
                energy.append([float(kdx * de) + min_e for kdx in range(int(n_jdx))])

        # add xes data to master xes list by interpolating to new energy range
        # first make new energy list
        if xes == []:
            for jdx, col in enumerate(x[0]):
                xes.append([0.0 for i in range(len(energy[jdx]))])
        # cycle through xes data
        for jdx in range(len(x[0])):
            # create temporary arrays for data
            e_temp = [e[i][jdx] for i in range(len(e))]
            x_temp = [x[i][jdx] for i in range(len(x))]
            # use 1d interpolation to map to new energy range
            f_interp = interp1d(e_temp, x_temp)
            # interpolate the new data
            #			print (min( e_temp ), max( e_temp ), min( energy[jdx] ), max( energy[jdx] ))
            xes_temp = f_interp(energy[jdx])

            # add the multiplicity-weighted XES to the master value
            for kdx in range(len(xes_temp)):
                xes[jdx][kdx] += mult_list[idx] * xes_temp[kdx] / mult_tot

    # open the summed file
    xes_data = open(sum_file, 'w')
    # write the header
    xes_data.write('K energy, K XS')
    if len(xes) > 2:
        xes_data.write(', L2,3 energy, L2,3 XS')
    if len(xes) > 4:
        xes_data.write(', M2,3 energy, M2,3 XS')
    xes_data.write('\n')
    # write in the data
    for jdx in range(len(energy[0])):
        for idx in range(len(energy)):
            if idx > 0:
                xes_data.write(', ')
            if len(energy[idx]) > jdx:
                xes_data.write(str(energy[idx][jdx]) + ', ' + str(xes[idx][jdx]))
            else:
                xes_data.write(', ')
        xes_data.write('\n')

    xes_data.close()
    return 0


####################################################################################################################################################
# MAIN SCRIPT ######################################################################################################################################
####################################################################################################################################################

# initialize directory list
calc_dir_list = []

# Rydberg to eV conversion
Ry_to_eV = 13.6056923

# initialize flags
no_sum = -1
no_nan = -1
no_sym = -1
no_tot = -1
no_xes = -1
updn = ''

# check that a file was specified
if len(sys.argv) > 1:
    for arg in sys.argv[1:]:
        # if the argument doesn't start with a dash, it should be the directory name
        if (arg[0] != '-'):
            # get rid of the trailing slash, if one was specified
            # this is just to make parsing the files in this directory simpler
            if arg[-1] == '/':
                calc_dir_list.append(arg[:-1])
            else:
                calc_dir_list.append(arg)
        # parse other arguments
        elif arg == '-All':
            # specify that we shouldn't sum the atomic types
            no_sum = 1
            print("Atoms of the same species will not be summed.")
        elif arg == '-NaN':
            # specify that we should average NaN values
            no_nan = 1
            print("NaN values will be excluded from clean DOS file.")
        elif arg == '-Sym':
            # specify that we want all the symmetries
            no_sym = 1
            print("All symmetries (s through f) will be calculated.")
        elif arg == '-Tot':
            # specify that we also want the total DOS
            no_tot = 1
            print("Total DOS will be calculated.")
        elif arg == '-up':
            # specify that we want the spin-up polarization
            updn = 'up'
            print("Using spin-up polarization.")
        elif arg == '-dn':
            # specify that we want the spin-down polarization
            updn = 'dn'
            print("Using spin-down polarization.")
        elif arg[:2] == '-X':
            # specify that we want to calculate the XES spectra
            no_xes = 1
            # see if we only want a particular symmetry
            if len(arg) > 2:
                if arg[2] == 'K':
                    no_xes = 2
                elif arg[2] == 'L':
                    no_xes = 3
                elif arg[2] == 'M':
                    no_xes = 4
            print("XES Spectra will be calculated.")

# parse the data in the directory
if calc_dir_list == []:
    print("Please specify a directory!")
    sys.exit()

# cycle through all directories given
for calc_dir in calc_dir_list:

    # let them know what is going on
    print("Processing DOS in directory", calc_dir)

    # initialize the lists for atomic information
    # atom_Z = []
    # atom_name = []
    # multiplicity = []
    # symmetry = []
    core_level = ()

    # get the name of this WIEN2k <case>
    case_name = calc_dir.split('/')[-1]
    # get the current directory (to return to later)
    start_dir = os.getcwd()

    # go to the <case> directory
    os.chdir(calc_dir)

    # get the multiplicity and symmetry of the atoms from the QTL file
    #	here the 'JATOM' flag identifies lines that contain the required information
    #	the multiplicity is the 4th entry, the symmetry is the last entry as a comma-separated string
    #	white space separates the main entries (the symmetry string has no spaces, just commas)
    return_val = get_params(case_name + '.qtl' + updn, 'JATOM', [3, -1], 0, delimiter=[' ', '=', '\n'])
    # make sure we have successfully read in the parameters from the QTL file
    if return_val == -1:
        sys.exit()
    # convert the multiplicity into an integer
    #print(return_val)
    multiplicity = [int(val) for val in return_val[0]]
    # leave the symmetry as a long string for easy handling
    symmetry = return_val[1]

    # get the atom names and atomic numbers from the STRUCT file
    #	here the 'NPT=' flag identifies lines that contain the required information
    #	the atom name is the 1st entry, the atomic number is the last
    return_val = get_params(case_name + '.struct', 'NPT=', [0, -1], idx_len=[3, 3])

    #print(return_val)
    # make sure we have successfully read in the parameters from the STRUCT file
    if return_val == -1:
        sys.exit()
    # here we don't need to convert to integers, since all we want for the atomic number is pattern matching
    atom_name = return_val[0]
    atom_Z = return_val[1]

    # get the number of core levels and binding energies from the SCFC file
    #	here the '#.ATOM'  flag identifies lines that contain the required information
    #	the number of core levels is the 3rd last entry
    # only do this if we are calculating XES
    if no_xes != -1:
        return_val = get_params(case_name + '.scfc' + updn, 'ATOM', [-3], flag_idx=0)

        # make sure we got something from the SCFC file
        if return_val == -1:
            print("Could not read in core levels from SCFC file!")
            sys.exit()

        # convert the number of core levels
        core_levels = ([int(val) for val in return_val[0]],)
        # initialize binding energies
        shell = [0.0 for val in return_val[0]]

        # look for K-shell binding energies
        if (no_xes == 1) or (no_xes == 2):
            return_val = get_params(case_name + '.scfc' + updn, '1S', [-2])
            # make sure we got something from the SCFC file
            if return_val == -1:
                print("Could not read in K-shell binding energies from SCFC file!")
                sys.exit()
            # make sure we add the binding energies to the proper atoms
            # I think every atom, even H, has a 1s core level (for H it is 0), but just in case...
            k_idx = 0
            for idx, val in enumerate(core_levels[0]):
                if val >= 1:
                    shell[idx] = Ry_to_eV * float(return_val[0][k_idx])
                    k_idx += 1
        # add K shell binding energies to tuple
        # we want to do this even if we aren't generating a K-shell XES
        core_levels += (shell,)
        # look for L-shell binding energies
        if (no_xes == 1) or (no_xes == 3):
            # get L2 energies
            return_val = get_params(case_name + '.scfc' + updn, '2P*', [-2])
            # make sure we got something from the SCFC file
            if return_val == -1:
                print("Could not read in L2-shell binding energies from SCFC file!")
                sys.exit()
            # initialize binding energies
            shell = [0.0 for val in range(len(core_levels[0]))]
            # make sure we add the binding energies to the proper atoms
            k_idx = 0
            for idx, val in enumerate(core_levels[0]):
                # for there to be a 2p* shell, we need to have at least 3 core levels (1s, 2s, 2p*)
                if val >= 3:
                    shell[idx] = Ry_to_eV * float(return_val[0][k_idx])
                    k_idx += 1
            # add L2 shell binding energies to tuple
            core_levels += (shell,)

            # get L3 energies
            return_val = get_params(case_name + '.scfc' + updn, '2P', [-2])
            # make sure we got something from the SCFC file
            if return_val == -1:
                print("Could not read in L3-shell binding energies from SCFC file!")
                sys.exit()
            # initialize binding energies
            shell = [0.0 for val in range(len(core_levels[0]))]
            # make sure we add the binding energies to the proper atoms
            k_idx = 0
            for idx, val in enumerate(core_levels[0]):
                # for there to be a 2p shell, we need to have at least 4 core levels (1s, 2s, 2p*, 2p)
                if val >= 4:
                    shell[idx] = Ry_to_eV * float(return_val[0][k_idx])
                    k_idx += 1
            # add L3 shell binding energies to tuple
            core_levels += (shell,)
        else:
            # add blanks for L2, and L3 binding energies
            shell = [0.0 for val in range(len(core_levels[0]))]
            core_levels += (shell, shell,)
        # look for M-shell binding energies
        if (no_xes == 1) or (no_xes == 4):
            # get M4 energies
            return_val = get_params(case_name + '.scfc' + updn, '3D*', [-2])
            # make sure we got something from the SCFC file
            if return_val == -1:
                print("Could not read in M4-shell binding energies from SCFC file!")
                sys.exit()
            # initialize binding energies
            shell = [0.0 for val in range(len(core_levels[0]))]
            # make sure we add the binding energies to the proper atoms
            k_idx = 0
            for idx, val in enumerate(core_levels[0]):
                # for there to be a 3d* shell, we need to have at least 8 core levels (1s, 2s, 2p*, 2p, 3s, 3p*, 3p, 3d*)
                if val >= 8:
                    shell[idx] = Ry_to_eV * float(return_val[0][k_idx])
                    k_idx += 1
            # add L3 shell binding energies to tuple
            core_levels += (shell,)

            # get M5 energies
            return_val = get_params(case_name + '.scfc' + updn, '3D', [-2])
            # make sure we got something from the SCFC file
            if return_val == -1:
                print("Could not read in M5-shell binding energies from SCFC file!")
                sys.exit()
            # initialize binding energies
            shell = [0.0 for val in range(len(core_levels[0]))]
            # make sure we add the binding energies to the proper atoms
            k_idx = 0
            for idx, val in enumerate(core_levels[0]):
                # for there to be a 3d* shell, we need to have at least 8 core levels (1s, 2s, 2p*, 2p, 3s, 3p*, 3p, 3d*, 3d)
                if val >= 9:
                    shell[idx] = Ry_to_eV * float(return_val[0][k_idx])
                    k_idx += 1
            # add L3 shell binding energies to tuple
            core_levels += (shell,)
        else:
            # add blanks for M4, and M5 binding energies
            shell = [0.0 for val in range(len(core_levels[0]))]
            core_levels += (shell, shell,)

    # check	to see that we have the correct data
    if (atom_Z == []) or (len(atom_Z) != len(multiplicity)) or (len(atom_Z) != len(symmetry)) or (
            len(atom_Z) != len(atom_name)):
        print("PROBLEM!!!")
        print("We have the following atom names : ", len(atom_name))
        print("We have the following atomic numbers : ", len(atom_Z))
        print("We have the following multiplicities : ", len(multiplicity))
        print("We have the following symmetries : ", len(symmetry))
        sys.exit()

    # cycle through atoms
    for atom in range(len(atom_Z)):
        # use TETRA to create DOS files
        tetra_flag = make_dos(case_name, str(atom + 1), symmetry[atom], no_sym)

        # if TETRA didn't work, stop the program
        if tetra_flag != 0:
            sys.exit()

        # clean up the DOS files
        clean_dos(case_name, atom_name[atom], no_nan)

        # FOR DEBUGGING ################
        # save the old DOS files
        # os_call = subprocess.call( ['mv', case_name + '.dos1ev', case_name + '.dos1ev_' + str( atom )] )
        # os_call = subprocess.call( ['mv', case_name + '.dos2ev', case_name + '.dos2ev_' + str( atom )] )
        ################################

        # scrap old TXSPEC files, if they exist
        for tx_stem in ['', '_K', '_L', '_M']:
            if os.path.exists(case_name + '.txspec' + updn + tx_stem) == True:
                rm_call = subprocess.call(['rm', case_name + '.txspec' + updn + tx_stem])
        # use XSPEC to make XES files, if necessary
        e1_k = None
        e1_l = None
        e1_m = None
        if (no_xes == 1) or (no_xes == 2):
            # check to see if k shell exists, if so generate it
            if core_levels[1][atom] < 0.0:
                e1_k = make_xes(case_name, str(atom + 1), 1)
        if (no_xes == 1) or (no_xes == 3):
            # check to see if l2 shell exists, if so generate it
            if core_levels[2][atom] < 0.0:
                e1_l = make_xes(case_name, str(atom + 1), 2)
        if (no_xes == 1) or (no_xes == 4):
            # check to see if d4 shell exists, if so generate it
            if core_levels[4][atom] < 0.0:
                e1_m = make_xes(case_name, str(atom + 1), 3)

        # clean XSPEC files, if necessary
        if no_xes > 0:
            binding_energy = [core_levels[i][atom] for i in range(1, 6)]
            clean_xes(case_name, atom_name[atom], 'XES', binding_energy, e1_k)

        ########
        # now do XAS
        # scrap old TXSPEC files, if they exist
        for tx_stem in ['', '_K', '_L', '_M']:
            if os.path.exists(case_name + '.txspec' + updn + tx_stem) == True:
                rm_call = subprocess.call(['rm', case_name + '.txspec' + updn + tx_stem])
        # use XSPEC to make XES files, if necessary
        if (no_xes == 1) or (no_xes == 2):
            # check to see if k shell exists, if so generate it
            if core_levels[1][atom] < 0.0:
                e1_k = make_xas(case_name, str(atom + 1), 1)
        if (no_xes == 1) or (no_xes == 3):
            # check to see if l2 shell exists, if so generate it
            if core_levels[2][atom] < 0.0:
                e1_l = make_xas(case_name, str(atom + 1), 2)
        if (no_xes == 1) or (no_xes == 4):
            # check to see if d4 shell exists, if so generate it
            if core_levels[4][atom] < 0.0:
                e1_m = make_xas(case_name, str(atom + 1), 3)

        # clean XSPEC files, if necessary
        if no_xes > 0:
            binding_energy = [core_levels[i][atom] for i in range(1, 6)]
            clean_xes(case_name, atom_name[atom], 'XAS', binding_energy, e1_k)

    # create total DOS if specified
    if no_tot == 1:
        # use TETRA to make the DOS files
        tetra_flag = make_dos(case_name, '0', 'tot', -1)
        clean_dos(case_name, 'TOT', no_nan)

    # check to see if we should sum the DOS for all equivalent atoms
    if no_sum == -1:
        # start a list of unique atoms
        # we could probably identify unique atoms by their atom_name, but using atom_Z is easier
        unique_atom_Z = []
        unique_atom_idx = []
        # cycle through the atoms
        for idx, Z in enumerate(atom_Z):
            # see if this atomic number exists in the list of unique numbers
            if unique_atom_Z.count(Z) == 0:
                # if so, add it to the list
                unique_atom_Z.append(Z)
                unique_atom_idx.append([idx])
            else:
                # otherwise, add this atom's index to the existing entry for this Z
                unique_atom_idx[unique_atom_Z.index(Z)].append(idx)

        # now gather the names of all the clean DOS files for each unique atom and add the data together
        for item in unique_atom_idx:
            # get the list of files to add together
            file_list = [case_name + '.DOS_' + atom_name[idx] for idx in item]
            # we want to write these files to a new file that just has the atom name, without any trailing identifiers
            sum_file = case_name + '.DOS_SUM_' + atom_name[item[0]].rstrip('0123456789')
            sum_mult = [multiplicity[idx] for idx in item]
            # add together the DOS files
            sum_dos(file_list, sum_file, sum_mult)

            # do this for XES as well, if applicable
            if no_xes > 0:
                file_list = [case_name + '.XES_' + atom_name[idx] for idx in item]
                sum_file = case_name + '.XES_SUM_' + atom_name[item[0]].rstrip('0123456789')
                # add together XES files, can use old multiplicity
                sum_xes(file_list, sum_file, sum_mult)

                file_list = [case_name + '.XAS_' + atom_name[idx] for idx in item]
                sum_file = case_name + '.XAS_SUM_' + atom_name[item[0]].rstrip('0123456789')
                # add together XES files, can use old multiplicity
                sum_xes(file_list, sum_file, sum_mult)

    # go back to starting directory
    os.chdir(start_dir)
