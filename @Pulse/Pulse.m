classdef Pulse
    %--- PULSE ------------------------------------------------------------
	%A housekeeping class to keep track of TRACK pulses and FITTED pulses,
	% as well as the bi-directional correspondence between the two sets.
	%
	% Matches tracked pulses (TRACK) to fitted pulses (FITTED), with a generalized
	% class MatchTrackFit to handle the two-way mapping between the two sets of
	% objects, which could contain the following categories:
	%	1) one-to-one matches
	%	2) missed objects from either set
	%	3) objects from one set mapped onto the same object in the other.
	%
	% Methods are implemented to handle manual adjustment to the correspondence
	% to resolve missed objects, as well as merged objects. To this end a graphical
	% output is also implemented under @pulse.graph.
	% 
	% Properties (private)
	%	fits - the complete array of all fits, does not have to exclude non-tracked
	%		cells/embryos
	%	fitsOI_ID - the fitID (unique ID code for each fitted pulse) of interest, i.e.
	%		those from cells with tracked pulses
	%	tracks - the array of TRACK pulses for this embryo
	%	tracks_mdf_file - the MDF filename
	%	cells - an array of CELLOBJ of data from cells found in this embryo
	%
	% Properties (public)
	% 	embryoID - the index of this embryo (see also LOAD_EDGE_SCRIPT)
	%	changes - a structure containing information about the manual changes
	%		performed on the PULSE object
	%
	% Methods
	% 
	% --- Construction methods ---
	%	Pulse - constructor object
	%	.match_pulse - performs matching between .tracks and .fits, and constructs
	%		the .map property from the resulting two-way mapping
	%	.categorize_mapping - creates/updates the .categories property from the .map
	%		includes the following fields:
	%		.one2one - one-to-one matches
	%		.merge - two Track per Fitted
	% 		.split - one Track per two Fitted
	%		.miss - Track with no Fitted
	% 		.add - Fitted with no Track
    %   .search_catID Search for the category_index (catID) of a pulse
    %       (track/fit).
    %   .cat - concatenate two Pulse objects (for example from different
    %       embryos)
	% --- Manual editing methods ---
    %   rename_embryoID - consistently renames the embryoID
	%	search_catID - given a tracked or fit from a PULSE object, find the index
	%		of that object within its current .cagetory
	%	removePulse - remove a track/fit and update the mapping
	%	createFitFromTrack - create an artificial FITTED from roughly where a TRACK
	%		was... uses MANUAL_FIT
	%	createTrackFromFit - create an artificial TRACK from roughly where a FITTED
	%		was...
	%	reassignFit - re-assign a FITTED to a TRACK, will only work if neither have
	% 		prior assignments
	%	read_changes - given a .change structure, edit current Pulse object
    %   adjust_centers - adjust the reference time used to construct
    %       dev_time
    %
    % --- Saving methods ---
    %   export_manual_fits - writes manually fitted parameters into a CSV
    %       file
	%	export_changes - writes all changes into a CSV
	% --- Display methods ---
	%	graph - Generates a 1x3 subplot of the TRACK/FIT/CELL
	% 	display - In-line display, reporting the number of objects and the quality of
	%		matching
    %
    % See also: CELLOBJ, FITTED, TRACK, FIND_ONE2ONE
	%
	% xies@mit.edu April 2013.

    properties
        
        fits 		% Complete set of FITTED pulses from ALL cells, including non-tracked and other embryos
        fitsOI_ID 	% fitID from only trakced cells in TRACKS
		fit_opt 	% The fitting option file for this embryo
        tracks 		% Set of TRACK pulses for this embryo
        tracks_mdf_file % The filename of the MDF file from which .tracks was loaded
        cells 		% An array of CELL of cells in this embryo (Contains the raw data)
        input      % input array with info about identifiy of each embryo
        next_fitID  % bookkeeping
        
        map			% The two-way mapping between Track and Fit
        match_thresh % The threshold of frames overlap above which a TRACK and a FITTED is matched (usually 1)
        categories %  A structure containing the different categories of matches
        
        embryoID
        changes
        
    end
    
    methods % Dynamic methods
        
% --------------------------- Constructor ---------------------------------

        function pulse = Pulse(tracks,filename,fits,opts,cells,input)
			%PULSE Constructor for the Pulse object (see main documentation)
			% Will not generate the .map property.
			%
			% USAGE: pulse = Pulse(tracks,mdf_filename,fits,fit_opt,cells,input);
            
            if nargin > 0 % empty case (for object arrays)
                
				% Make sure that tracks are from the same embryo
				if any( [tracks.embryoID] ~= tracks(1).embryoID ),
					error('Error making Pulse object: Tracks need to come from the same embryo');
				end
                
                if strcmp(class(tracks),'Track')
                    pulse.fits = fits;
                    pulse.tracks = tracks;
                else
                    pulse.fits = tracks;
                    pulse.tracks = fits;
                end
                
                fitsOI_ID = [fits( ... filter out non-tracked cells
                    ismember( [fits.stackID], [tracks.stackID] )).fitID];
                
                pulse.tracks_mdf_file = filename;
                pulse.input = input;
                pulse.fitsOI_ID = fitsOI_ID;
                pulse.next_fitID = fits.get_fitID(fitsOI_ID(1)).embryoID*100000;
                pulse.fit_opt = opts( fits.get_fitID(fitsOI_ID(1)).embryoID );

				% Take only cells from the same embryoID
                cells = cells( [cells.embryoID] == tracks(1).embryoID );
                pulse.cells = cells;
                
                save([fileparts(pulse.tracks_mdf_file) '/pulse_raw.mat'], 'pulse');
                
            end % non-empty constructor
            
        end % Constructor
        
% ------------------------ Fitted handling --------------------------------
        
        assign_datafield(pulse,name);
        align_fits(pulse,name,measurement);
        interpolate_traces(pulse,name,dt);
        retrace(pulse, opts);
        
% ------------------------ Pulse measurements -----------------------------
        
        M = get_corrected_measurement(pulse,meas,input);
        myosin_persistence = get_myosin_persistence(fits);
        [perc,varargout] = percent_overlap(fits,cells);
        
% -------------------- Neighbor-neighbor measurements ---------------------
        
        fits = find_near_fits(pulse,neighbor_def);
        f = find_non_edge(fits,cells);
        num_near = get_num_near(pulse,neighbor_definition,window);
        
        fits_bs = bootstrap_cluster_label(fits);
        [fits_bs,cells_bs] = simulate_pulsing(fits,cells,freqHat);
        
% ---------------------- Cell-level analysis ------------------------------

        % Pulsing analysis
        first_fits = get_first_fit(pulse);
        [freq,center] = get_frequency(pulse);
        [freq,neighbor_count] = estimate_pulsing_params(pulse);
        [adj,nodes] = get_pulsing_trajectories(pulse);
        [adj,nodes] = get_pulse_transition_graph(pulse);
        W = get_pulse_transition_matrix(pulse);
        
%         [fits_bs,cells_bs] = monte_carlo_stackID(pulse)
        
% --------------------------- Array handling ------------------------------
        
        pulse = cat(pulse1,pulse2);
        
% --------------------------- Mapping -------------------------------------
        
        pulse = categorize_mapping(pulse);
        pulse = match_pulse(pulse,threshold);
        catID = search_catID(pulse,type,pulseID);
        
% --------------------- Edit pulse/tracks ---------------------------------
        
        pulse = removePulse(pulse,type,pulseID);
        pulse = createTrackFromFit(pulse,fitID);
        pulse = createFitFromTrack(pulse,trackID,opt);
        pulse = reassignFit(pulse,fitID,newTrackID);
		pulse = read_changes( pulse, changes );
        
% ----------------------- Find object for robust exporting ----------------
        
        obj = find_nearest_object(pulse,obj_type,cx,cy,ct);
        
% ----------------------- Saving / exporting ------------------------------

        export_manual_fits(pulse);
        [cx,cy,ct] = get_xyt(pulse);
        
        function export_changes( pulse )
            %EXPORT_CHANGES
            % Export all .changes to a .mat file
            changes = pulse.changes;
            save( [fileparts(pulse.tracks_mdf_file), '/', 'changes.mat'], 'changes');
            % export manual changes
            pulse.export_manual_fits;
        end % export_changes
        
% ---------------------- graph/display ------------------------------------
        
        varargout = graph(pulse,cat,ID,axes_handle);
        binary = make_binary_sequence(pulse);
%         disp(pulse)
        
% ---------------------- Edit embryo-level parameters ---------------------
        
        function pulse = adjust_dev_time(pulse,input)
            % ADJUST_DEV_TIME
            % 
            % Adjust the .center and .dev_time of CellObj and Fitted
            % objects with respect to new reference development itme.
            
            if numel(pulse) > 1
                error('Only one embryo please.')
            end

            old_tref = find(pulse.cells(1).dev_time == 0);
            
            pulse.input = input;
            new_tref = input.tref;
            dt = input.dt;

            pulse.cells = pulse.cells.adjust_dev_time(old_tref,new_tref,dt);
            pulse.fits = pulse.fits.adjust_dev_time(old_tref,new_tref,dt);
            
        end % Adjust adjust_dev_time
        
        function pulse = rename_embryoID(pulse,embryoID)
            % Rename all the FITTED and CELLOBJ from an old embryoID into a
            % new embryoID.
            %
            % USAGE: pulse = pulse.rename_embryoID(newID)
            
            %@TODO:  ALSO UPDATE the MAP!!!!!
            
            old_embryoID = pulse.embryoID;
            pulse.embryoID = embryoID;
            pulse.tracks = pulse.tracks.rename_embryoID(embryoID);
            pulse.fits = pulse.fits.rename_embryoID(embryoID);
            pulse.cells = pulse.cells.rename_embryoID(embryoID);
            
            for i = 1:numel(pulse.fitsOI_ID)
                
                fID = pulse.fitsOI_ID(i);
                base = 10.^floor(log10(fID) - log10(old_embryoID));
                fID = fID - old_embryoID*base + embryoID*base;
                pulse.fitsOI_ID(i) = fID;
                
            end
            
        end % rename_embryoID
        
    end % Dynamic methods
    
end
