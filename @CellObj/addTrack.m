function cellobj = addTrack(cellobj,trackID)
%@CellObj.addTrack Add a trackID to a cell
cellobj.trackID = [ [cellobj.trackID] trackID];
cellobj.num_tracks = cellobj.num_tracks + 1;
end
