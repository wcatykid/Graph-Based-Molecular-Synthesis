/*
 *  This file is part of esynth.
 *
 *  esynth is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  esynth is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with esynth.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _TIMED_LIKE_VALUE_CONTAINER_GUARD
#define _TIMED_LIKE_VALUE_CONTAINER_GUARD 1


#include <string>
#include <algorithm>
#include <vector>
#include <sstream>


class TimedLikeValueContainer
{
  private:
    struct TimedStruct
    {
        unsigned time;
        std::string data;
        
        TimedStruct(std::string s, unsigned t) : data(s), time(t) {}
        
        inline bool operator==(const TimedStruct& rhs)
        {
            return this->data == rhs.data;
        }

        inline bool equals(const std::string& that) const
        {
            return this->data == that;
        }
    };

    // The container for the records
    std::vector<TimedStruct>* table;
    static unsigned time_rel_UPPERBOUND;
    static unsigned time_ABS_UPPERBOUND;

  public:
    inline int size() const { return table->size(); }
    inline bool empty() const { return table->empty(); }

    static void SetThresholds(unsigned relUB, unsigned absUB)
    {
        time_rel_UPPERBOUND = relUB;
        time_ABS_UPPERBOUND = absUB;
    }

    TimedLikeValueContainer()
    {
        table = new std::vector<TimedStruct>();
    }
   
    // Kill all entries
    ~TimedLikeValueContainer()
    {
        table->clear();
        delete table;
    }
	
    bool contains(const std::string& that) const
    {
        for (std::vector<TimedStruct>::const_iterator it = table->begin(); it != table->end(); it++)
        {
            if (it->equals(that)) return true;
        }
        
        return false;
    }
	
    // We assume a call to contains has already been made (to limit redundant calls)
    void add(const std::string& that, unsigned theTime)
    {
        table->push_back(TimedStruct(that, theTime));
    }
	
    // We assume a call to contains has already been made (to limit redundant calls)
    void PurgeAndRenumber()
    {
        std::vector<TimedStruct>* newtable = new std::vector<TimedStruct>();

        // Copy the elements from the old list to the new list purging 
        for (std::vector<TimedStruct>::const_iterator it = table->begin(); it != table->end(); it++)
        {
            // If the given data is 'old'; that is, it's time is > than the relative upperbound
            if (it->time <= time_rel_UPPERBOUND)
            {
                newtable->push_back(TimedStruct(it->data, it->time + (time_ABS_UPPERBOUND - time_rel_UPPERBOUND)));
            }
        }

        //
        // The new table will replace the old.
        //
        table->clear();
        delete table;
        table = newtable;
    }

    std::string toString() const
    {
        std::ostringstream oss;
        
        for (std::vector<TimedStruct>::const_iterator it = table->begin(); it != table->end(); it++)
        {
            oss << it->data << "(" << it->time << ") ";
        }
        
        return oss.str();
    }
    
    friend std::ostream& operator<<(std::ostream& os, const TimedLikeValueContainer& container)
    {
        os << container.toString();
        
        return os;
    }
};

#endif
