function pulse = categorize_mapping(pulse)
%CATEGORIZE_MAPPING Quantify the different types of matches between
% FITTED pulses and TRACK pulses.
%
% USAGE: pulse = categorize_mapping(pulse)
% Updates the following class properties:
%	pulse.categories.one2one - pulses with one FITTED and one TRACK
%	pulse.categories.merge - one TRACK with multiple FITTED
%	pulse.categories.split - multiple TRACK with single FITTED
%	pulse.categories.miss - TRACK with no FITTED
%	pulse.categories.add - FITTED with no TRACK
%
% xies@mit.edu Feb 2013.

nbm = pulse.map;
fit = pulse.fits.get_fitID(pulse.fitsOI_ID);
track = pulse.tracks;

% -------- Quantify merges --
[trackID,fitID] = find_merges(nbm.dictTrackFit);
[ matches.merge( 1:numel(trackID) ).trackID ] = deal(trackID{:});
[ matches.merge( 1:numel(trackID) ).fitID ] = deal(fitID{:});
% Annotate fit/track with merges
for i = 1:numel(trackID)
    [ fit( ismember([fit.fitID], fitID{i}) ).category ] = deal('merge');
    [ track( ismember([track.trackID], trackID{i}) ).category ] = deal('merge');
end

% ------- Quantify splits --
[fitID,trackID] = find_merges(nbm.dictFitTrack);
[ matches.split( 1:numel(trackID) ).fitID ] = deal(fitID{:});
[ matches.split( 1:numel(trackID) ).trackID ] = deal(trackID{:});
% Annotate fit/track with splits
for i = 1:numel(trackID)
    [ fit( ismember([fit.fitID], fitID{i}) ).category ] = deal('split');
    [ track( ismember([track.trackID], trackID{i}) ).category ] = deal('split');
end

% ------- Quantify missed/added --
matchedTF_trackID = nbm.dictTrackFit.keys; %fitID
matchedFT_fitID = nbm.dictFitTrack.keys; %trackID
% misses
trackID = ...
    [ track(~ismember([track.trackID],cell2mat(matchedTF_trackID))).trackID ];
% Annotate track with misses
matches.miss.trackID = [];
matches.miss.fitID = [];
for i = 1:numel(trackID)
    matches.miss(i).trackID = trackID(i);
    matches.miss(i).fitID = [];
    track( [track.trackID] == trackID(i) ).category = 'miss';
end

% adds
matches.add.trackID = [];
matches.add.fitID = [];
fitID = ...
    [ fit(~ismember([fit.fitID],cell2mat(matchedFT_fitID))).fitID ];
for i = 1:numel(fitID)
    matches.add(i).fitID = fitID(i);
    matches.add(i).trackID = [];
    fit( [fit.fitID] == fitID(i) ).category = 'add';
end

% -- Quantify one-to-one matches --
% Find all tracks within the fit-fitted cells, and not belonging to
% merge/split/miss
one2one_origins = [track( ...
    ~ismember([track.trackID], unique([matches.merge.trackID matches.split.trackID])) ...
    & ~ismember([track.trackID], [matches.miss.trackID]) ...
    & ismember([track.trackID], cell2mat(nbm.dictTrackFit.keys)) ).trackID];
% initialize
[one2one(1:numel(one2one_origins)).trackID] = deal([]);
[one2one(1:numel(one2one_origins)).fitID] = deal([]);
% Assign
for i = 1:numel(one2one_origins)
    one2one(i).trackID = one2one_origins(i);
    one2one(i).fitID = nbm.dictTrackFit( one2one_origins(i) );
end
% Annotate one2one matches ontp fit/track structures
[ track( ismember([track.trackID],[one2one.trackID]) ).category] = deal('one2one');
[ fit( ismember([fit.fitID],[one2one.fitID]) ).category] = deal('one2one');
% Collect into structure
matches.one2one = one2one;

% Delete empty fields
matches = delete_empty(matches);

pulse.fits = fit;
pulse.tracks = track;
pulse.categories = matches;

% -- Subfunctions of categorize_mapping -- %
    function match = delete_empty(match)
        if isempty( [match.one2one.trackID] )
            match = rmfield(match,'one2one');
        end
        if isempty( [match.merge.trackID] )
            match = rmfield(match,'merge');
        end
        if isempty( [match.split.trackID] )
            match = rmfield(match,'split');
        end
        if isempty( [match.miss.trackID] )
            match = rmfield(match,'miss');
        end
        if isempty( [match.add.fitID] )
            match = rmfield(match,'add');
        end
    end
% --- End categorize_mapping subfunctions
end % Categorize_mapping
