#ifndef MI_EVENT_HANDLER_MACROS
#define MI_EVENT_HANDLER_MACROS

// re-implement wx-like macros for our versions

#undef BEGIN_EVENT_TABLE
#undef END_EVENT_TABLE
#undef EVT_MENU
#undef EVT_UPDATE_UI

#define BEGIN_EVENT_TABLE(klas, parent) { \
        MIEventHandler *_dispatcher = klas; \
        MIChildEventHandlerFtor *ftor = 0;

#define END_EVENT_TABLE() }

#define EVT_MENU(id, handler) \
    _dispatcher->registerActionHandler(id, # handler, ftor);

#define EVT_UPDATE_UI(id, handler) \
    _dispatcher->registerUpdateHandler(id, # handler, ftor);

#endif // MI_EVENT_HANDLER_MACROS
