%for par in parameters:
   ${par["fortrantype"]}${par["extratype"]}${par["len"]} :: ${par["name"]}${par["dimension"]} ${par["value"]}
%endfor
## vim: filetype=mako
