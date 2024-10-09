<%
numvars = len(variables)
maxnamelen = max(len(variable["name"]) for variable in variables)
def rank(variable):
    return len(variable["shape"])
maxrank = max(rank(variable) for variable in variables)
%>
        integer, parameter :: numvars    =  ${numvars}
        integer, parameter :: maxnamelen =  ${maxnamelen}
        integer, parameter :: maxrank    =  ${maxrank}
%for var in variables:
   character(len=maxnamelen),parameter,public ::  mnem_${var['name'].ljust(maxnamelen)} = '${var['name'].ljust(maxnamelen)}'
%endfor
        character(len=maxnamelen),dimension(numvars),parameter :: mnemonics= (/ &
%for i, var in enumerate(variables):
         mnem_${var['name'].ljust(maxnamelen)}${"," if i != len(variables)-1 else ''} & ! ${i+1}
%endfor
        /)

## vim: filetype=mako 
