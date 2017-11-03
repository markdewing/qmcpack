
// Temporary converter from PackedBuffer to OOMPI_Message
// Can be removed when OOMPI is removed.

#include <Utilities/PackedBuffer.h>
#include <OOMPI/Message.h>

inline OOMPI_Message toMessage(PackedBuffer &p)
{
  return OOMPI_Message(p.buf, p.bufSize);
}
