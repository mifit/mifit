#ifndef MI_CONFIGBASE_H
#define MI_CONFIGBASE_H

#include <string>
#ifdef _WIN32
#undef GetProfileInt
#undef GetProfileString
#undef WriteProfileInt
#undef WriteProfileString
#endif

class MIConfig {
  public:
    static MIConfig *Instance() {
      return _instance;
    }
    
    virtual ~MIConfig() {
    }
    
    // note: for better or worse, our implementation of config has
    // two overloads of read for each type:
    //    one which does not take a default argument, and which does not change the value if it's not in the file
    //    one which does
    // we avoid reimplementation by introducing the assignDefaultValue variable
    virtual bool Read(const std::string& key, std::string& value, const std::string& defaultValue, bool assignDefaultValue = true) = 0;
    virtual bool Read(const std::string& key, long* value, long defaultValue, bool assignDefaultValue = true) = 0;
    virtual bool Read(const std::string& key, bool* value, bool defaultValue, bool assignDefaultValue = true) = 0;
    virtual bool Read(const std::string& key, double* value, double defaultValue, bool assignDefaultValue = true) = 0;


    virtual bool Write(const std::string& key, const std::string& value) = 0;
    virtual bool Write(const std::string& key, long value) = 0;
    virtual bool Write(const std::string& key, bool value) = 0;
    virtual bool Write(const std::string& key, double value) = 0;

    virtual void Flush() = 0;

    bool Read(const std::string& key, std::string& value);
    bool Read(const std::string& key, long* value);
    bool Read(const std::string& key, bool* value);
    bool Read(const std::string& key, double* value);

    int GetProfileInt(const std::string& group, const std::string& key, int defaultValue = int());
    std::string GetProfileString(const std::string& group, const std::string& key, const std::string& defaultValue = "");

    void WriteProfileInt(const std::string& group, const std::string& key, int value);
    void WriteProfileString(const std::string& group, const std::string& key, const std::string& value);

    virtual bool HasKey(const std::string &key) = 0;

  protected:
    static MIConfig* _instance;
    MIConfig();
};

#endif
