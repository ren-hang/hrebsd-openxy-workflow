function Phase = ValidatePhase(Phase)
%Validate Material Detection
MaterialsList = GetMaterialsList(2);
if ~all(ismember(Phase,MaterialsList))
    invalidMats = unique(Phase(~ismember(Phase,MaterialsList)));
    % LOCAL PATCH (headless): map unknown phase names (e.g. 'notindexed') to one named material,
    % the same choice protocol/Brady make by hand in the 'Material detection failed' dialog.
    hm = getenv('OPENXY_HEADLESS_MATERIAL');
    if ~isempty(hm) && ismember(hm,MaterialsList)
        warning('ValidatePhase:headless','Mapping phases {%s} -> ''%s'' (OPENXY_HEADLESS_MATERIAL)',strjoin(invalidMats,', '),hm);
        Phase(~ismember(Phase,MaterialsList)) = {hm};
        return;
    end
    er = errordlg(['Auto material detection failed. "' strjoin(invalidMats,', ') '" not found in list of known materials'],'Material Detection');
    uiwait(er,5)
    for i = 1:length(invalidMats)
        op = questdlg(['Select an option for "' invalidMats{i} '":'],'Material Detection Failed','Select Existing Material','Create a New Material','Cancel','Select Existing Material');
        inds = strcmp(Phase,invalidMats{i});
        while true
            switch op
                case 'Select Existing Material'
                    Materials = GetMaterialsList(3);
                    [index, ok] = listdlg('PromptString','Select a Material','ListString',Materials,'SelectionMode','single','Name','Material Selection');
                    if ok
                        Phase(inds) = {Materials{index}};
                        break;
                    else
                        op = 'Cancel';
                    end
                case 'Create a New Material'
                    material = NewMaterialGUI;
                    if material ~= 0
                        Phase(inds) = {material};
                        break;
                    else
                        op = 'Cancel';
                    end
                case 'Cancel'
                    er = warndlg('Material selection failed. Select a new Scan File.','Material Selection');
                    uiwait(er)
                    Phase = {};
                    break;
            end
        end
    end
end
