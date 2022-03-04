function str_ver = check_MTEX_version(current_version)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
str_ver = getMTEXpref('version');
required = split(extractAfter(str_ver, ' '),'.');
current = split(extractAfter(current_version,'-'),'.');
r_major = str2double(required{1});
r_minor = str2double(required{2});
c_major = str2double(current{1});
c_minor = str2double(current{2});

assert(r_major == c_major, "Your MTEX version" + current + " is incompatable with" +required);
assert(r_minor >=c_minor, "Your MTEX version" + current + " is incompatable with" + required);

end

