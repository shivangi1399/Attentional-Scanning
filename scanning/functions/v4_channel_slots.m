function [rows, slots, nBad] = v4_channel_slots(labels, animal, nChTot)
% Shared helper -- a VERBATIM extraction of the local subfunction of the same
% name in erp_latency_wave.m, so that ongoing_excitability_wave.m,
% wave_behaviour_coupling.m and anything else reading clean_lfp.mat /
% clean_mua.mat uses the SAME channel mapping. Getting it wrong rotates the
% array by one electrode and silently averages channel 43 onto 44, so it now
% exists on the path rather than only inside one script.
%
% NOTE ON THE DUPLICATE. erp_latency_wave.m keeps its own local copy and MATLAB
% resolves locals before the path, so that script is unaffected by this file and
% continues to use its own. The two are identical today. If the upstream label
% convention ever changes, BOTH must be updated -- or the local removed so this
% one takes over. Same arrangement as functions/elec_rf_deg.m.
%
% Map this session's V4 channels onto the canonical 1..nChTot slots BY LABEL
% NUMBER, not by sorted position: sessions do not all contain all 64 V4
% channels (hermes 60/64, klecks 62/64) and which ones are missing varies, so
% position would put a different electrode in row i in different sessions. This
% mirrors how the phase pipeline builds its 64-channel array
% (Phase_analysis/masters_code/Phase_combine_sessions.m): each present channel
% goes to its own slot, absent channels are left NaN, never compacted. Slot i
% here is channel i in phase_progression.mat.
%
%   hermes   V4-n -> slot n        (labels 1..64,   blank ->  64)
%   klecks   V4-n -> slot n - 64   (labels 65..128, blank -> 128)
%
% The klecks offset is -64, not -63: clean_lfp.mat already holds the incremented
% numbering, not mapping_lfp.m's pre-increment 64..127. Verified against
% phase_progression.mat on three exact matches -- hermes zero-trial channels
% {48, 64} equal the slots absent from every session; klecks V4-109 is absent
% from all 25 sessions -> slot 45, the only klecks channel with zero trials;
% klecks V4-128 appears in one session -> slot 64, which has 68 trials against
% ~1400 elsewhere. With -63 every klecks channel would sit one slot off.
%
%   rows   row indices into D.label / D.trial for the channels kept
%   slots  canonical slot for each of those rows
%   nBad   V4 labels numbered outside 1..nChTot (should be 0)

labels = labels(:);
isV4 = startsWith(strtrim(labels), 'V4-');
V4i  = find(isV4);
nums = str2double(erase(strtrim(labels(isV4)), 'V4-'));   % ' 1' -> 1, '' -> NaN
switch lower(animal)
    case 'hermes'
        nums(isnan(nums)) = 64;
        slots = nums;
    case 'klecks'
        nums(isnan(nums)) = 128;
        slots = nums - 64;
    otherwise
        error('Unknown animal: %s', animal);
end
keep  = slots >= 1 & slots <= nChTot & isfinite(slots);
nBad  = sum(~keep);
rows  = V4i(keep);
slots = slots(keep);
if numel(unique(slots)) ~= numel(slots)
    error('Duplicate canonical slots in %s -- label convention is wrong.', animal);
end
end
