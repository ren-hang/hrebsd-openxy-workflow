function varargout = uigetdir(varargin)
% Headless shim: any GUI dialog becomes a logged error instead of blocking forever.
msg = ''; if nargin>0, a = varargin{1}; if iscell(a), msg = strjoin(a(:)',' | '); elseif ischar(a), msg = a; end; end
if nargin>1 && ischar(varargin{2}), msg = [msg ' || ' varargin{2}]; end
error('Headless:uigetdir','Dialog uigetdir() blocked in headless run: %s', msg);
end

