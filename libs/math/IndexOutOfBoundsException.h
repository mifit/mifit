#ifndef mi_math_IndexOutOfBoundsException_h
#define mi_math_IndexOutOfBoundsException_h

namespace mi
{
    namespace math
    {

        class IndexOutOfBoundsException
        {
        private:
            char *detail;
        public:
            IndexOutOfBoundsException();
            IndexOutOfBoundsException(const char *detail);
            IndexOutOfBoundsException(const IndexOutOfBoundsException &ex);
            IndexOutOfBoundsException&operator=(const IndexOutOfBoundsException &ex);
            virtual ~IndexOutOfBoundsException();
            virtual const char *what() const;
        };

    }
}


#endif // ifndef mi_math_IndexOutOfBoundsException_h
