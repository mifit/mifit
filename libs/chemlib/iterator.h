#ifndef MI_ITERATOR_H
#define MI_ITERATOR_H

//Base class for iterators
template<class T>
class MIIterBase
{
public:
    MIIterBase()
    {
    }

    virtual MIIterBase<T> *Clone() const = 0;
    virtual ~MIIterBase()
    {
    }

    virtual MIIterBase<T>&First() = 0;
    virtual MIIterBase<T>&Last() = 0;
    virtual MIIterBase<T>&Start() = 0;

    virtual MIIterBase<T>&operator++() = 0;
    virtual MIIterBase<T>&operator--() = 0;

    virtual T&operator*() const = 0;
    virtual T*operator->() const = 0;
    virtual operator bool() const
    {
        return false;
    }

    virtual operator T*() const
    {
        return operator->();
    }

};


// wrapper class around MIIterBase
//   handles automatic dereferencing of iterators
//   takes ownership of iterator and automatically deletes it

template<class T>
class MIIter
{
public:
    MIIter(MIIterBase<T> *rhs) : _i(rhs)
    {
    }

    MIIter(const MIIter<T> &rhs) : _i(rhs.Clone())
    {
    }

    MIIter<T>&operator=(const MIIter<T> &rhs)
    {
        if (this != &rhs)
        {
            delete _i;
            _i = rhs.Clone();
        }
        return *this;
    }

    MIIter&operator=(const MIIter<T> *rhs)
    {
        if (this != rhs)
        {
            delete _i;
            _i = rhs->Clone();
        }
        return *this;
    }

    MIIter&operator=(const MIIterBase<T> &rhs)
    {
        if (this != &rhs)
        {
            delete _i;
            _i = rhs.Clone();
        }
        return *this;
    }

    //Must return an MIIterBase*, rather than a MIIter*, otherwise we can't assign it
    MIIterBase<T> *Clone() const
    {
        return _i->Clone();
    }

    ~MIIter()
    {
        delete _i;
    }

    MIIter<T>&First()
    {
        _i->First();
        return *this;
    }

    MIIter&Last()
    {
        _i->Last();
        return *this;
    }

    MIIter&Start()
    {
        _i->Start();
        return *this;
    }

    MIIter&operator++()
    {
        _i->operator++();
        return *this;
    }

    MIIter&operator--()
    {
        _i->operator--();
        return *this;
    }

    T&operator*() const
    {
        return _i->operator*();
    }

    T*operator->() const
    {
        return _i->operator->();
    }

    // automatic casting operators
    operator T*() const {
        return operator->();
    }

    operator T& () const {
        return *(operator->());
    }

    operator bool() const {
        return _i->operator bool();
    }

private:
    MIIter() : _i(0)
    {
    }

    MIIterBase<T> *_i;
};

template<class T>
class MIDoublyLinkedListIter : public MIIterBase<T>
{
public:
    MIDoublyLinkedListIter(T *start) : _start(start), _current(start)
    {
    }

    virtual MIIterBase<T> *Clone() const
    {
        MIDoublyLinkedListIter<T> *i;
        if (_current != NULL)
        {
            i = new MIDoublyLinkedListIter<T>(_current);
        }
        else
        {
            i = new MIDoublyLinkedListIter<T>(_start);
        }
        return i;
    }

    virtual MIIterBase<T>&First()
    {
        if (_start != NULL)
        {
            if (_current == NULL)
            {
                _current = _start;
            }
            T *t = _current->prev();
            while (t != NULL)
            {
                _current = t;
                t = _current->prev();
            }
        }
        return *this;
    }

    virtual MIIterBase<T>&Last()
    {
        if (_start != NULL)
        {
            if (_current == NULL)
            {
                _current = _start;
            }
            T *t = _current->next();
            while (t != NULL)
            {
                _current = t;
                t = _current->next();
            }
        }
        return *this;
    }

    virtual MIIterBase<T>&Start()
    {
        _current = _start;
        return *this;
    }

    virtual MIIterBase<T>&operator++()
    {
        if (_current != NULL)
        {
            _current = _current->next();
        }
        return *this;
    }

    virtual MIIterBase<T>&operator--()
    {
        if (_current != NULL)
        {
            _current = _current->prev();
        }
        return *this;
    }

    virtual T&operator*() const
    {
        return *_current;
    }

    virtual T*operator->() const
    {
        return _current;
    }

    virtual operator bool() const
    {
        return _start != NULL && _current != NULL;
    }

private:
    T *_start;
    T *_current;
};

#endif // ifndef MI_ITERATOR_H

