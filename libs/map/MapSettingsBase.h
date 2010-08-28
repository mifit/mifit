#ifndef mifit_map_MapSettingsBase_h
#define mifit_map_MapSettingsBase_h

class MapSettingsBase
{
public:
    /**
     * Countering method to use - see CONTOUR_METHOD_TYPE
     */
    int ContourMethod;
    /**
     * Radius for the blob contouring method.
     */
    float BlobRadius;
    /**
     * Radius of the contouring box.
     */
    float Radius;
    /**
     * Map levels to contour.
     */
    float MapLevel[5];
    /**
     * If a map level is used or not
     */
    bool MapLevelOn[5];
    /**
     * Color of each map level.
     */
    int MapLevelColor[5];
    /**
     * Width of map contouring lines in pixels
     */
    float maplinewidth;

    int m_radiusmax;

    float mapmin;
    float mapmax;

    MapSettingsBase();
    virtual ~MapSettingsBase()
    {
    }

};
#endif // ifndef mifit_map_MapSettingsBase_h
