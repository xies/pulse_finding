function cellobj = removeTrack(cellobj,trackID)
%@CellObj.addTrack Add a trackID from a cell
cellobj.trackID([cellobj.trackID] == trackID) = [];
cellobj.num_tracks = cellobj.num_tracks - 1;
end