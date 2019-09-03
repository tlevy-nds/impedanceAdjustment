classdef MessageUpdater < handle
    properties
        messagestr
    end
    
    methods
        function obj = MessageUpdater()
            obj.messagestr = '';
        end
        
        function update_message(obj, newmessage)        
            backspaces = repmat('\b', size(obj.messagestr));
            obj.messagestr = newmessage;
            fprintf([backspaces obj.messagestr]);
        end
        
        function clear(obj)
            obj.messagestr = '';
        end
    end
end
