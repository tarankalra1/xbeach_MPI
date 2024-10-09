  !  DO NOT EDIT THIS FILE
  !  But edit variables.f90 and scripts/generate.py
  !  Compiling and running is taken care of by the Makefile

  <%
formats = ["a", "afafafaf", "afaiaaa", "aaafaf", "aafaf", "afaaa", "aafa", "aaaf",
"afafa", "afaf", "afa", "aaf", "iiiii", "aiaiafa",
"aiaiaf", "aiaiaia", "aiaiai", "aiafaf", "aiafa", "aaaiai",
"aaiai", "aiai", "aiaaa", "aaai",
"aii", "aai", "aaaa", "aaa", "ai", "aa",  "illll", "af",
"aia", "ia", "fa","aaia","aiaa", "aiaia"]
%>
  !
  ! Options for destiantion in writelog
  ! 's' = screen
  ! 'l' = log file
  ! 'e' = error file
  ! 'w' = warning file
  !
  ! Combinations also allowed, f.i.
  ! 'le' = log file and error file
  ! 'el' ditto
  ! 'sel' = screen, log file and error file
  !
  interface writelog
%for format in formats:
     module procedure writelog_${format}
%endfor
  end interface writelog

## vim: filetype=mako

