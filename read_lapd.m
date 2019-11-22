function read_lapd()

file_name = "/Volumes/DATA/LAPD/RF_April22/100_BaPSF_check_plasma_profile.hdf5";

info = h5info(file_name);

% data = info.Groups(2).Groups(3).Datasets(5);

% ch6_name = '/Raw data + config/SIS crate/Mach_antennacurrent_longertime [Slot 7: SIS 3302 ch 6]';
% ch6_name = '/Raw data + config/SIS crate/Bdot_only_antenna [Slot 5: SIS 3302 ch 1]';
% ch6_name = '/Raw data + config/SIS crate/Isat [Slot 5: SIS 3302 ch 1]';


% disp('other channel names are ...');
% 
% for i=1:numel(data)
%     
%     disp(data(i).Name);
%     
% end

% exract some subsection from the data set ...
%
% info.Groups(2).Groups(3).Datasets(1).Dataspace(1)
% info.Groups(2).Groups(3).Datasets(5).Dataspace
% tells us that the data size is
% [120832 x 10670]
% I don't know which dimension is which, but we can extract pieces like ...

startX = 20000;
startY = 1;
nX = 1;
nY = 4268;

ch6_data = h5read(file_name,ch6_name, [startX startY], [nX nY]);

% to get all the data for that channel simply use ...

ch6_data = h5read(file_name,ch6_name);



end
  