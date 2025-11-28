% First get the template file (F_co2dist) read in

ncid = netcdf.open('F_co2dist.nc');

varid = netcdf.inqVarID(ncid,'F_co2dist');

data = netcdf.getVar(ncid,0);
netcdf.close(ncid);



plast = zeros([100 100 150]);
year1 = 1:150;

% For co2dist data
% Use estimate of all plastic produced from Geyer et al. 2017
% Compound growth in production of 8.4%/y

mp1 = 2. * 1.084.^(year1);
year = [year1 year1(end)];
mp = [mp1 mp1(end)];

% convert MT to particles using mass conversion from Eriksen et al. (2014)
% 355.4E8 gMP = 485E10 particles
  % And per year to per second, and per m2
mp = mp.*1e6.*485e10./355.4e8./365./24./60./60./0.356684e15;

% Now distribute spatially according to F_co2dist
for x = 1:1:151
    plast(:,:,x) = mp(x).*data(:,:,11);
end

ncid = netcdf.create('newfile.nc', 'NOCLOBBER');
dimid1 = netcdf.defDim(ncid, 'time', 151);
dimid2 = netcdf.defDim(ncid, 'latitude', 100);
dimid3 = netcdf.defDim(ncid, 'longitude', 100);
ncdf_time = netcdf.defVar(ncid, 'time', 'float', dimid1);
ncdf_mp = netcdf.defVar(ncid, 'O_mp','float',[dimid2 dimid3 dimid1]);
varid = netcdf.getConstant('Global');

netcdf.endDef(ncid);

netcdf.putVar(ncid,ncdf_time, year);
netcdf.putVar(ncid,ncdf_mp, plast);

netcdf.reDef(ncid);

netcdf.putAtt(ncid, varid , 'title', 'MP flux'); 
netcdf.putAtt(ncid, varid, 'history', 'Created by Karin Kvale'); %Update

netcdf.putAtt(ncid, ncdf_time, 'standard_name', 'time');
netcdf.putAtt(ncid, ncdf_time, 'long_name', 'time');
netcdf.putAtt(ncid, ncdf_time, 'axis', 'T');
netcdf.putAtt(ncid, ncdf_time , 'units', 'year');
netcdf.putAtt(ncid, ncdf_time, 'missing_value', 9.969209968386869e+36);
netcdf.putAtt(ncid, ncdf_time , 'fill_value', 9.969209968386869e+36);

netcdf.putAtt(ncid, ncdf_mp , 'standard_name', 'O_MP');
netcdf.putAtt(ncid, ncdf_mp , 'long_name', 'microplastic');
netcdf.putAtt(ncid, ncdf_mp , 'units', 'particles per m^-2 s^-1');
netcdf.putAtt(ncid, ncdf_mp , 'missing_value', 9.969209968386869e+36);
netcdf.putAtt(ncid, ncdf_mp , 'fill_value', 9.969209968386869e+36);

netcdf.close(ncid);


